/******************************************************************************
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
******************************************************************************/

#ifndef Minisat_SolverTypes_h
#define Minisat_SolverTypes_h

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace Minisat {

//=============================================================================
// Variables

typedef int Var;
constexpr Var var_Undef(-1);

struct VariableOrder {
  std::vector<float> &activity;

  bool operator()(Var x, Var y) const { return activity[x] > activity[y]; }
};

//=============================================================================
// Literals

class Lit {
protected:
  int x;

  explicit constexpr Lit(int x, int b) : x(x ^ b) {}

public:
  Lit() = default;
  constexpr Lit(const Lit &) = default;
  explicit constexpr Lit(Var var, bool sign = false)
      : x(var + var + (int)sign) {}

  constexpr bool operator==(Lit p) const { return x == p.x; }
  constexpr bool operator!=(Lit p) const { return x != p.x; }
  constexpr bool operator<(Lit p) const {
    return x < p.x; // '<' makes p, ~p adjacent in the ordering.
  }

  constexpr Lit operator~() const { return Lit(x, 1); }
  Lit operator^(bool b) const { return Lit(x, (int)b); }
  constexpr bool sign() const { return x & 1; }
  constexpr Var var() const { return x >> 1; }
  constexpr int toInt() const { return x; }
};

constexpr Lit lit_Undef(var_Undef, false);
constexpr Lit lit_Error(var_Undef, true);

//=============================================================================
// Lifted booleans

class lbool {
  uint8_t value;

public:
  constexpr explicit lbool(uint8_t v) : value(v) {}

  constexpr lbool() : value(0) {}
  explicit constexpr lbool(bool x) : value(!x) {}

  constexpr bool operator==(lbool b) const {
    return (value & b.value & 2) != 0 || value == b.value;
  }

  constexpr bool operator!=(lbool b) const { return !(*this == b); }
  constexpr lbool operator^(bool b) const {
    return lbool((uint8_t)(value ^ (uint8_t)b));
  }

  constexpr bool isTrue() const { return value == 0; }
  constexpr bool isFalse() const { return value == 1; }
};

constexpr lbool l_True((uint8_t)0);
constexpr lbool l_False((uint8_t)1);
constexpr lbool l_Undef((uint8_t)2); // 3 is also undef

//============================================================================
// Clause

class Clause;
typedef Clause *CRef;

class Clause {
private:
  std::vector<Lit> lits;
  float act;

public:
  static constexpr CRef UNDEF = NULL;

  Clause(std::vector<Lit> &lits, bool learnt)
      : lits(std::move(lits)), act(learnt ? 0.0f : -1.0f) {}

  constexpr int size() const { return lits.size(); }
  constexpr bool learnt() const { return act >= 0.0f; }
  constexpr const std::vector<Lit> &literals() const { return lits; }

  Lit &operator[](int i) { return lits[i]; }
  Lit operator[](int i) const { return lits[i]; }

  float &activity() {
    assert(act >= 0.0f);
    return act;
  }
};

struct ClauseSatisfied {
  std::vector<lbool> &assigns;

  inline bool operator()(CRef c) const {
    for (Lit p : c->literals())
      if ((assigns[p.var()] ^ p.sign()) == l_True)
        return true;
    return false;
  }
};

//============================================================================
// Watcher

struct Watcher {
  CRef cref;
  Lit blocker;

  constexpr Watcher(CRef cr, Lit p) : cref(cr), blocker(p) {}

  constexpr bool operator==(const Watcher &w) const { return cref == w.cref; }
  constexpr bool operator!=(const Watcher &w) const { return cref != w.cref; }
};

//============================================================================
// VarData

struct VarData {
  CRef reason;
  int level;
};

} // namespace Minisat

#endif
