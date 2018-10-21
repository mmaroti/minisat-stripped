/****************************************************************************************[Solver.h]
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

#ifndef Minisat_Solver_h
#define Minisat_Solver_h

#include <vector>

#include "minisat/mtl/Heap.h"
#include "SolverTypes.h"


namespace Minisat {

class Solver {
public:
  // Constructor:
  Solver();

  // Problem specification:
  Lit addLiteral(bool polarity = true, bool decision = true);
  bool giveClause(std::vector<Lit> &ps);      // empties the source vector
  bool addClause();                           // make the solver contraditory
  bool addClause(Lit p);                      // add unit clause
  bool addClause(Lit p, Lit q);               // add binary clause
  bool addClause(Lit p, Lit q, Lit r);        // add ternary clause
  bool addClause(const std::vector<Lit> &ps); // ret false if contradictory

  // Solving:
  bool simplify();                 // removes already satisfied clauses
  bool solve();                    // search without assumptions
  bool solve(Lit p);               // use a single assumption
  bool solve(Lit p, Lit q);        // use two assumptions
  bool solve(Lit p, Lit q, Lit r); // use tree assumptions
  bool solve(const std::vector<Lit> &assumps);
  bool okay() const; // check if solver not in contradictory
  void interrupt();  // ask solver to return cleanly
  void clearInterrupt();

  // Variable mode:
  //
  void    setPolarity    (Var v, bool b); // Declare which polarity the decision heuristic should use for a variable. Requires mode 'polarity_user'.
  void    setDecisionVar (Var v, bool b); // Declare if a variable should be eligible for selection in the decision heuristic.

  // Read state:
  //
  lbool   value      (Var x) const;       // The current value of a variable.
  lbool   value      (Lit p) const;       // The current value of a literal.
  lbool   modelValue (Lit p) const;       // The value of a literal in the last model. The last call to solve must have been satisfiable.
  int     nAssigns   ()      const;       // The current number of assigned literals.
  int     nClauses   ()      const;       // The current number of original clauses.
  int     nLearnts   ()      const;       // The current number of learnt clauses.
  int     nVars      ()      const;       // The current number of variables.
  int     nFreeVars  ()      const;

  // Extra results: (read-only member variable)
  //
  std::vector<lbool> model;     // If problem is satisfiable, this vector contains the model (if any).
  std::vector<Lit>   conflict;  // If problem is unsatisfiable (possibly under assumptions),
                                // this vector represent the final conflict clause expressed in the assumptions.

  // Mode of operation:
  //
  double    var_decay;
  double    clause_decay;
  double    random_var_freq;
  double    random_seed;
  bool      luby_restart;
  bool      rnd_pol;            // Use random polarities for branching heuristics.
  bool      rnd_init_act;       // Initialize variable activities with a small random value.
  double    garbage_frac;       // The fraction of wasted memory allowed before a garbage collection is triggered.

  int       restart_first;      // The initial restart limit.                                                                (default 100)
  double    restart_inc;        // The factor with which the restart limit is multiplied in each restart.                    (default 1.5)
  double    learntsize_factor;  // The intitial limit for learnt clauses is a factor of the original clauses.                (default 1 / 3)
  double    learntsize_inc;     // The limit for learnt clauses is multiplied with this factor each restart.                 (default 1.1)

  int       learntsize_adjust_start_confl;
  double    learntsize_adjust_inc;

  // Statistics: (read-only member variable)
  //
  long solves, starts, decisions, rnd_decisions, propagations, conflicts;
  long dec_vars, clauses_literals, learnts_literals, max_literals, tot_literals;

protected:
  // Helper structures:
  //
  struct VarData { CRef reason; int level; };
  static inline VarData mkVarData(CRef cr, int l){ VarData d = {cr, l}; return d; }

  struct VarOrderLt {
      const std::vector<float>&  activity;
      bool operator () (Var x, Var y) const { return activity[x] > activity[y]; }
      VarOrderLt(const std::vector<float>&  act) : activity(act) { }
  };

  // Solver state:
  //
  bool                ok;               // If FALSE, the constraints are already unsatisfiable. No part of the solver state may be used!
  std::vector<CRef>   clauses;          // List of problem clauses.
  std::vector<CRef>   learnts;          // List of learnt clauses.
  double              cla_inc;          // Amount to bump next clause with.
  std::vector<float>  activity;         // A heuristic measurement of the activity of a variable.
  double              var_inc;          // Amount to bump next variable with.
  std::vector<std::vector<Watcher>> watches; // 'watches[lit]' is a list of constraints watching 'lit' (will go there if literal becomes true).
  std::vector<lbool>  assigns;          // The current assignments.
  std::vector<bool>   polarity;         // The preferred polarity of each variable.
  std::vector<char>   decision;         // Declares if a variable is eligible for selection in the decision heuristic.
  std::vector<Lit>    trail;            // Assignment stack; stores all assigments made in the order they were made.
  std::vector<int>    trail_lim;        // Separator indices for different decision levels in 'trail'.
  std::vector<VarData> vardata;         // Stores reason and level for each variable.
  int                 qhead;            // Head of queue (as index into the trail -- no more explicit propagation queue in MiniSat).
  int                 simpDB_assigns;   // Number of top-level assignments since last execution of 'simplify()'.
  long                simpDB_props;     // Remaining number of propagations that must be made before next execution of 'simplify()'.
  std::vector<Lit>    assumptions;      // Current set of assumptions provided to solve by the user.
  Heap<VarOrderLt>    order_heap;       // A priority queue of variables ordered with respect to the variable activity.

  // Temporaries (to reduce allocation overhead)
  std::vector<bool> analyze_seen;
  std::vector<Lit> analyze_stack;
  std::vector<Lit> analyze_toclear;
  std::vector<Lit> addclause_temp;

  double              max_learnts;
  double              learntsize_adjust_confl;
  int                 learntsize_adjust_cnt;

  // Resource contraints:
  //
  bool                asynch_interrupt;

  // Main internal methods:
  //
  std::vector<Watcher>& occurences(Lit p) { return watches[p.toInt()]; }             // occurence list lookup
  void     insertVarOrder   (Var x);                                                 // Insert a variable in the decision order priority queue.
  Lit      pickBranchLit    ();                                                      // Return the next decision variable.
  void     newDecisionLevel ();                                                      // Begins a new decision level.
  void     uncheckedEnqueue (Lit p, CRef from = Clause::UNDEF);                      // Enqueue a literal. Assumes value of literal is undefined.
  bool     enqueue          (Lit p, CRef from = Clause::UNDEF);                      // Test if fact 'p' contradicts current state, enqueue otherwise.
  CRef     propagate        ();                                                      // Perform unit propagation. Returns possibly conflicting clause.
  void     cancelUntil      (int level);                                             // Backtrack until a certain level.
  void     analyze          (CRef confl, std::vector<Lit>& out_learnt, int& out_btlevel); // (bt = backtrack)
  void     analyzeFinal     (Lit p, std::vector<Lit>& out_conflict);                 // COULD THIS BE IMPLEMENTED BY THE ORDINARIY "analyze" BY SOME REASONABLE GENERALIZATION?
  bool     litRedundant     (Lit p, unsigned int abstract_levels);                            // (helper method for 'analyze()')
  lbool    search           (int nof_conflicts);                                     // Search for a given number of conflicts.
  lbool    solve_           ();                                                      // Main solve method (assumptions given in 'assumptions').
  void     reduceDB         ();                                                      // Reduce the set of learnt clauses.
  void     rebuildOrderHeap ();

  // Maintaining Variable/Clause activity:
  void varDecayActivity() { var_inc *= (1.0 / var_decay); }
  void varBumpActivity(Var v);
  void claDecayActivity() { cla_inc *= (1.0 / clause_decay); }
  void claBumpActivity(Clause &c);

  // Operations on clauses:
  //
  void     attachClause     (CRef cr);               // Attach a clause to watcher lists.
  void     detachClause     (CRef cr); // Detach a clause to watcher lists.
  void     removeClause     (CRef cr);               // Detach and free a clause.
  bool     locked           (const Clause& c) const; // Returns TRUE if a clause is a reason for some implication in the current state.

  // Misc:
  //
  int      decisionLevel    ()      const; // Gives the current decisionlevel.
  unsigned int abstractLevel(Var x) const; // Used to represent an abstraction of sets of decision levels.
  CRef     reason           (Var x) const { return vardata[x].reason; }
  int      level            (Var x) const { return vardata[x].level; }

  // Static helpers:
  //

  // Returns a random float 0 <= x < 1. Seed must never be 0.
  inline double drand() {
    random_seed *= 1389796;
    int q = (int)(random_seed / 2147483647);
    random_seed -= (double)q * 2147483647;
    return random_seed / 2147483647;
  }

  // Returns a random integer 0 <= x < size. Seed must never be 0.
  inline int irand(int size) { return (int)(drand() * size); }
};

// Problem specification:

inline bool Solver::addClause() {
  addclause_temp.clear();
  return giveClause(addclause_temp);
}

inline bool Solver::addClause(Lit p) {
  addclause_temp.clear();
  addclause_temp.push_back(p);
  return giveClause(addclause_temp);
}

inline bool Solver::addClause(Lit p, Lit q) {
  addclause_temp.clear();
  addclause_temp.push_back(p);
  addclause_temp.push_back(q);
  return giveClause(addclause_temp);
}

inline bool Solver::addClause(Lit p, Lit q, Lit r) {
  addclause_temp.clear();
  addclause_temp.push_back(p);
  addclause_temp.push_back(q);
  addclause_temp.push_back(r);
  return giveClause(addclause_temp);
}

inline bool Solver::addClause(const std::vector<Lit> &ps) {
  addclause_temp = ps;
  return giveClause(addclause_temp);
}

// Solving:

inline bool Solver::solve() {
  assumptions.clear();
  return solve_() == l_True;
}

inline bool Solver::solve(Lit p) {
  assumptions.clear();
  assumptions.push_back(p);
  return solve_() == l_True;
}

inline bool Solver::solve(Lit p, Lit q) {
  assumptions.clear();
  assumptions.push_back(p);
  assumptions.push_back(q);
  return solve_() == l_True;
}

inline bool Solver::solve(Lit p, Lit q, Lit r) {
  assumptions.clear();
  assumptions.push_back(p);
  assumptions.push_back(q);
  assumptions.push_back(r);
  return solve_() == l_True;
}

inline bool Solver::solve(const std::vector<Lit> &assumps) {
  assumptions = assumps;
  return solve_() == l_True;
}

inline bool Solver::okay() const { return ok; }

inline void Solver::interrupt() { asynch_interrupt = true; }

inline void Solver::clearInterrupt() { asynch_interrupt = false; }

// Activity

inline void Solver::insertVarOrder(Var x) {
  if (!order_heap.inHeap(x) && decision[x])
    order_heap.insert(x);
}

inline void Solver::varBumpActivity(Var v) {
  if ((activity[v] += var_inc) > 1e20f) {
    for (int i = 0; i < nVars(); i++)
      activity[i] *= 1e-20f;
    var_inc *= 1e-20;
  }

  if (order_heap.inHeap(v))
    order_heap.decrease(v);
}

inline void Solver::claBumpActivity(Clause &c) {
  if ((c.activity() += cla_inc) > 1e20f) {
    for (Clause *learnt : learnts)
      learnt->activity() *= 1e-20f;
    cla_inc *= 1e-20;
  }
}

// NOTE: enqueue does not set the ok flag! (only public methods do)
inline bool     Solver::enqueue         (Lit p, CRef from)      { return value(p) != l_Undef ? value(p) != l_False : (uncheckedEnqueue(p, from), true); }
inline bool     Solver::locked          (const Clause& c) const { return value(c[0]) == l_True && reason(c[0].var()) != Clause::UNDEF && reason(c[0].var()) == &c; }
inline void     Solver::newDecisionLevel()                      { trail_lim.push_back(trail.size()); }

inline int      Solver::decisionLevel ()      const   { return trail_lim.size(); }
inline unsigned int Solver::abstractLevel (Var x) const { return 1 << (level(x) & (8 * sizeof(unsigned int) - 1)); }
inline lbool    Solver::value         (Var x) const   { return assigns[x]; }
inline lbool    Solver::value         (Lit p) const   { return assigns[p.var()] ^ p.sign(); }
inline lbool    Solver::modelValue    (Lit p) const   { return model[p.var()] ^ p.sign(); }
inline int      Solver::nAssigns      ()      const   { return trail.size(); }
inline int      Solver::nClauses      ()      const   { return clauses.size(); }
inline int      Solver::nLearnts      ()      const   { return learnts.size(); }
inline int      Solver::nVars         ()      const   { return vardata.size(); }
inline int      Solver::nFreeVars     ()      const   { return (int)dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]); }
inline void     Solver::setPolarity   (Var v, bool b) { polarity[v] = b; }
inline void     Solver::setDecisionVar(Var v, bool b)
{
    if      ( b && !decision[v]) dec_vars++;
    else if (!b &&  decision[v]) dec_vars--;

    decision[v] = b;
    insertVarOrder(v);
}

}

#endif
