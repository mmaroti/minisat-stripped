/***************************************************************************************[Solver.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <math.h>
#include <algorithm>

#include "Solver.h"

using namespace Minisat;

// Constructor:

Solver::Solver()
    : // User parameters
      var_decay(0.95), clause_decay(0.999), random_var_freq(0),
      random_seed(91648253), luby_restart(true),
      restart_first(100), restart_inc(2),

      // Internal parameters
      learntsize_factor(1.0 / 3.0), learntsize_inc(1.1),
      learntsize_adjust_start_confl(100), learntsize_adjust_inc(1.5),

      // Statistics
      solves(0), starts(0), decisions(0), rnd_decisions(0), propagations(0),
      conflicts(0), clauses_literals(0), learnts_literals(0),
      max_literals(0), tot_literals(0),

      // State
      ok(true), cla_inc(1), var_inc(1), watches(), qhead(0), simpDB_assigns(-1),
      simpDB_props(0), order_heap({activity}), asynch_interrupt(false) {}

// Problem specification:

Lit Solver::addLiteral() {
  int v = nVars();
  int l = std::max(Lit(v, false).toInt() + 1, Lit(v, true).toInt() + 1);
  if (watches.size() < l)
    watches.resize(l);

  assigns.push_back(l_Undef);
  vardata.push_back({Clause::UNDEF, 0});
  activity.push_back(0.0);
  analyze_seen.push_back(false);
  insertVarOrder(v);

  return Lit(v, true);
}

bool Solver::takeClause(std::vector<Lit> &ps) {
  assert(decisionLevel() == 0);
  if (!ok)
    return false;

  // Check if clause is satisfied and remove false/duplicate literals:
  std::sort(ps.begin(), ps.end());
  Lit p = lit_Undef;
  auto i = ps.begin();
  auto j = ps.begin();
  auto end = ps.end();
  while (i != end) {
    if (value(*i) == l_True || *i == ~p) {
      return true;
    } else if (value(*i) != l_False && *i != p) {
      *j = p = *i;
      ++j;
    }
    ++i;
  }
  ps.erase(j, ps.end());

  if (ps.empty())
    return ok = false;
  else if (ps.size() == 1) {
    uncheckedEnqueue(ps[0]);
    return ok = (propagate() == Clause::UNDEF);
  }

  CRef cr = new Clause(ps, false);
  clauses.push_back(cr);
  attachClause(cr);
  return true;
}

void Solver::attachClause(CRef cr) {
    const Clause& c = *cr;
    assert(c.size() > 1);
    occurences(~c[0]).push_back(Watcher(cr, c[1]));
    occurences(~c[1]).push_back(Watcher(cr, c[0]));
    if (c.learnt()) learnts_literals += c.size();
    else            clauses_literals += c.size(); }


void Solver::detachClause(CRef cr) {
    const Clause& c = *cr;
    assert(c.size() > 1);

    // remove(occurences(~c[0]), Watcher(cr, c[1]));
    std::vector<Watcher> &occ0 = occurences(~c[0]);
    *std::find(occ0.begin(), occ0.end(), Watcher(cr, c[1])) = occ0.back();
    occ0.pop_back();

    // remove(occurences(~c[1]), Watcher(cr, c[0]));
    std::vector<Watcher> &occ1 = occurences(~c[1]);
    *std::find(occ1.begin(), occ1.end(), Watcher(cr, c[0])) = occ1.back();
    occ1.pop_back();

    if (c.learnt()) learnts_literals -= c.size();
    else            clauses_literals -= c.size(); }


void Solver::removeClause(CRef cr) {
    Clause& c = *cr;
    detachClause(cr);
    // Don't leave pointers to free'd memory!
    if (locked(c)) vardata[c[0].var()].reason = Clause::UNDEF;
    delete cr;
}

// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int level) {
    if (decisionLevel() > level){
        for (int c = trail.size()-1; c >= trail_lim[level]; c--){
            Var      x  = trail[c].var();
            assigns [x] = l_Undef;
            insertVarOrder(x); }
        qhead = trail_lim[level];
        trail.resize(trail_lim[level]);
        trail_lim.resize(level);
    } }


//=================================================================================================
// Major methods:


Lit Solver::pickBranchLit()
{
    Var next = var_Undef;

    // Random decision:
    if (drand() < random_var_freq && !order_heap.empty()){
        next = order_heap[irand(order_heap.size())];
        if (value(next) == l_Undef)
            rnd_decisions++; }

    // Activity based decision:
    while (next == var_Undef || value(next) != l_Undef)
        if (order_heap.empty()){
            next = var_Undef;
            break;
        }else
            next = order_heap.removeMin();

    return next == var_Undef ? lit_Undef : Lit(next, drand() < 0.5);
}


/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : std::vector<Lit>&) (out_btlevel : int&)  ->  [void]
|
|  Description:
|    Analyze conflict and produce a reason clause.
|
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the
|        rest of literals. There may be others from the same level though.
|
|________________________________________________________________________________________________@*/
void Solver::analyze(CRef confl, std::vector<Lit>& out_learnt, int& out_btlevel)
{
    int pathC = 0;
    Lit p     = lit_Undef;

    // Generate conflict clause:
    //
    out_learnt.push_back(lit_Undef);      // (leave room for the asserting literal)
    int index   = trail.size() - 1;

    do{
        assert(confl != Clause::UNDEF); // (otherwise should be UIP)
        Clause& c = *confl;

        if (c.learnt())
            claBumpActivity(c);

        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];

            if (!analyze_seen[q.var()] && level(q.var()) > 0){
                varBumpActivity(q.var());
                analyze_seen[q.var()] = true;
                if (level(q.var()) >= decisionLevel())
                    pathC++;
                else
                    out_learnt.push_back(q);
            }
        }

        // Select next clause to look at:
        while (!analyze_seen[trail[index--].var()]);
        p     = trail[index+1];
        confl = reason(p.var());
        analyze_seen[p.var()] = false;
        pathC--;

    }while (pathC > 0);
    out_learnt[0] = ~p;

    // Simplify conflict clause:
    //
    int i, j;
    analyze_toclear = out_learnt;

    unsigned int abstract_level = 0;
    for (i = 1; i < out_learnt.size(); i++)
        abstract_level |= abstractLevel(out_learnt[i].var()); // (maintain an abstraction of levels involved in conflict)

    for (i = j = 1; i < out_learnt.size(); i++)
        if (reason(out_learnt[i].var()) == Clause::UNDEF || !litRedundant(out_learnt[i], abstract_level))
            out_learnt[j++] = out_learnt[i];

    max_literals += out_learnt.size();
    out_learnt.resize(j);
    tot_literals += out_learnt.size();

    // Find correct backtrack level:
    //
    if (out_learnt.size() == 1)
        out_btlevel = 0;
    else{
        int max_i = 1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 2; i < out_learnt.size(); i++)
            if (level(out_learnt[i].var()) > level(out_learnt[max_i].var()))
                max_i = i;
        // Swap-in this literal at index 1:
        Lit p             = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
        out_btlevel       = level(p.var());
    }

    for (auto const& elem : analyze_toclear) {
        analyze_seen[elem.var()] = false;
    }
}


// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, unsigned int abstract_levels)
{
    analyze_stack.clear(); analyze_stack.push_back(p);
    int top = analyze_toclear.size();
    while (analyze_stack.size() > 0){
        assert(reason(analyze_stack.back().var()) != Clause::UNDEF);
        Clause& c = *reason(analyze_stack.back().var()); analyze_stack.pop_back();

        for (int i = 1; i < c.size(); i++){
            Lit p  = c[i];
            if (!analyze_seen[p.var()] && level(p.var()) > 0){
                if (reason(p.var()) != Clause::UNDEF && (abstractLevel(p.var()) & abstract_levels) != 0){
                    analyze_seen[p.var()] = true;
                    analyze_stack.push_back(p);
                    analyze_toclear.push_back(p);
                }else{
                    for (int j = top; j < analyze_toclear.size(); j++)
                        analyze_seen[analyze_toclear[j].var()] = false;
                    analyze_toclear.resize(top);
                    return false;
                }
            }
        }
    }

    return true;
}


/*_________________________________________________________________________________________________
|
|  analyzeFinal : (p : Lit)  ->  [void]
|
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
|    stores the result in 'out_conflict'.
|________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Lit p, std::vector<Lit>& out_conflict)
{
    out_conflict.clear();
    out_conflict.push_back(p);

    if (decisionLevel() == 0)
        return;

    analyze_seen[p.var()] = true;

    for (int i = trail.size()-1; i >= trail_lim[0]; i--){
        Var x = trail[i].var();
        if (analyze_seen[x]){
            if (reason(x) == Clause::UNDEF){
                assert(level(x) > 0);
                out_conflict.push_back(~trail[i]);
            }else{
                Clause& c = *reason(x);
                for (int j = 1; j < c.size(); j++)
                    if (level(c[j].var()) > 0)
                        analyze_seen[c[j].var()] = true;
            }
            analyze_seen[x] = false;
        }
    }

    analyze_seen[p.var()] = false;
}


void Solver::uncheckedEnqueue(Lit p, CRef from)
{
    assert(value(p) == l_Undef);
    assigns[p.var()] = lbool(!p.sign());
    vardata[p.var()] = {from, decisionLevel()};
    trail.push_back(p);
}


/*_________________________________________________________________________________________________
|
|  propagate : [void]  ->  [Clause*]
|
|  Description:
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
|    otherwise Clause::UNDEF.
|
|    Post-conditions:
|      * the propagation queue is empty, even if there was a conflict.
|________________________________________________________________________________________________@*/
CRef Solver::propagate()
{
    CRef    confl     = Clause::UNDEF;
    int     num_props = 0;

    while (qhead < trail.size()){
        Lit            p   = trail[qhead++];     // 'p' is enqueued fact to propagate.
        std::vector<Watcher>&  ws  = occurences(p);
        std::vector<Watcher>::iterator i, j, end;
        num_props++;

        for (i = j = ws.begin(), end = ws.end();  i != end;){
            // Try to avoid inspecting the clause:
            Lit blocker = i->blocker;
            if (value(blocker) == l_True){
                *j++ = *i++; continue; }

            // Make sure the false literal is data[1]:
            CRef     cr        = i->cref;
            Clause&  c         = *cr;
            Lit      false_lit = ~p;
            if (c[0] == false_lit)
                c[0] = c[1], c[1] = false_lit;
            assert(c[1] == false_lit);
            i++;

            // If 0th watch is true, then clause is already satisfied.
            Lit     first = c[0];
            Watcher w     = Watcher(cr, first);
            if (first != blocker && value(first) == l_True){
                *j++ = w; continue; }

            // Look for new watch:
            for (int k = 2; k < c.size(); k++)
                if (value(c[k]) != l_False){
                    c[1] = c[k]; c[k] = false_lit;
                    occurences(~c[1]).push_back(w);
                    goto NextClause; }

            // Did not find watch -- clause is unit under assignment:
            *j++ = w;
            if (value(first) == l_False){
                confl = cr;
                qhead = trail.size();
                // Copy the remaining watches:
                while (i < end)
                    *j++ = *i++;
            }else
                uncheckedEnqueue(first, cr);

        NextClause:;
        }
        ws.erase(j, end);
    }
    propagations += num_props;
    simpDB_props -= num_props;

    return confl;
}


/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
namespace {
    struct reduceDB_lt {
        bool operator () (CRef x, CRef y) {
            return x->size() > 2 && (y->size() == 2 || x->activity() < y->activity());
        }
    };
}

void Solver::reduceDB()
{
    int     i, j;
    float  extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity

    std::sort(learnts.begin(), learnts.end(), reduceDB_lt());
    // Don't delete binary or locked clauses. From the rest, delete clauses from the first half
    // and clauses with activity smaller than 'extra_lim':
    for (i = j = 0; i < learnts.size(); i++){
        Clause& c = *learnts[i];
        if (c.size() > 2 && !locked(c) && (i < learnts.size() / 2 || c.activity() < extra_lim))
            removeClause(learnts[i]);
        else
            learnts[j++] = learnts[i];
    }
    learnts.resize(j);
}

void Solver::rebuildOrderHeap()
{
    std::vector<Var> vs;
    for (Var v = 0; v < nVars(); v++)
        if (value(v) == l_Undef)
            vs.push_back(v);
    order_heap.build(vs);
}

/*_________________________________________________________________________________________________
|
|  simplify : [void]  ->  [bool]
|
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
bool Solver::simplify()
{
    assert(decisionLevel() == 0);

    if (!ok || propagate() != Clause::UNDEF)
        return ok = false;

    if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
        return true;

    // Remove satisfied clauses:
    ClauseSatisfied pred = {assigns};
    learnts.erase(std::remove_if(learnts.begin(), learnts.end(), pred),
                  learnts.end());
    clauses.erase(std::remove_if(clauses.begin(), clauses.end(), pred),
                  clauses.end());

    simpDB_assigns = nAssigns();
    simpDB_props   = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

    return true;
}


/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
|
|  Description:
|    Search for a model the specified number of conflicts.
|    NOTE! Use negative value for 'nof_conflicts' indicate infinity.
|
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool Solver::search(int nof_conflicts)
{
    assert(ok);
    int         backtrack_level;
    int         conflictC = 0;
    std::vector<Lit>    learnt_clause;
    starts++;

    for (;;){
        CRef confl = propagate();
        if (confl != Clause::UNDEF){
            // CONFLICT
            conflicts++; conflictC++;
            if (decisionLevel() == 0) return l_False;

            learnt_clause.clear();
            analyze(confl, learnt_clause, backtrack_level);
            cancelUntil(backtrack_level);

            if (learnt_clause.size() == 1){
                uncheckedEnqueue(learnt_clause[0]);
            }else{
                Lit p = learnt_clause[0];
                CRef cr = new Clause(learnt_clause, true);
                learnts.push_back(cr);
                attachClause(cr);
                claBumpActivity(*cr);
                uncheckedEnqueue(p, cr);
            }

            varDecayActivity();
            claDecayActivity();

            if (--learntsize_adjust_cnt == 0){
                learntsize_adjust_confl *= learntsize_adjust_inc;
                learntsize_adjust_cnt    = (int)learntsize_adjust_confl;
                max_learnts             *= learntsize_inc;
            }

        }else{
            // NO CONFLICT
            if (nof_conflicts >= 0 && (conflictC >= nof_conflicts || asynch_interrupt)){
                // Reached bound on number of conflicts:
                cancelUntil(0);
                return l_Undef; }

            // Simplify the set of problem clauses:
            if (decisionLevel() == 0 && !simplify())
                return l_False;

            if (learnts.size()-nAssigns() >= max_learnts)
                // Reduce the set of learnt clauses:
                reduceDB();

            // New variable decision:
            decisions++;
            Lit next = pickBranchLit();

            if (next == lit_Undef)
                // Model found:
                return l_True;

            // Increase decision level and enqueue 'next'
            newDecisionLevel();
            uncheckedEnqueue(next);
        }
    }
}


/*
  Finite subsequences of the Luby-sequence:

  0: 1
  1: 1 1 2
  2: 1 1 2 1 1 2 4
  3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
  ...


 */

static double luby(double y, int x){

    // Find the finite subsequence that contains index 'x', and the
    // size of that subsequence:
    int size, seq;
    for (size = 1, seq = 0; size < x+1; seq++, size = 2*size+1);

    while (size-1 != x){
        size = (size-1)>>1;
        seq--;
        x = x % size;
    }

    return pow(y, seq);
}

lbool Solver::solve_()
{
    model.clear();
    if (!ok) return l_False;

    solves++;

    max_learnts               = nClauses() * learntsize_factor;
    learntsize_adjust_confl   = learntsize_adjust_start_confl;
    learntsize_adjust_cnt     = (int)learntsize_adjust_confl;
    lbool   status            = l_Undef;

    // Search:
    int curr_restarts = 0;
    while (status == l_Undef){
        double rest_base = luby_restart ? luby(restart_inc, curr_restarts) : pow(restart_inc, curr_restarts);
        status = search(rest_base * restart_first);
        if (asynch_interrupt) break;
        curr_restarts++;
    }

    if (status == l_True){
        // Extend & copy model:
        model.resize(nVars());
        for (int i = 0; i < nVars(); i++) model[i] = value(i);
    }else if (status == l_False)
        ok = false;

    cancelUntil(0);
    return status;
}
