/***********************************************************************************[SolverTypes.h]
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


#ifndef Minisat_SolverTypes_h
#define Minisat_SolverTypes_h

#include <cassert>
#include <algorithm>
#include <vector>

#include "minisat/mtl/Vec.h"

namespace Minisat {

//=================================================================================================
// Variables, literals, lifted booleans, clauses:


// NOTE! Variables are just integers. No abstraction here. They should be chosen from 0..N,
// so that they can be used as array indices.

typedef int Var;
#define var_Undef (-1)

class Lit {
protected:
  int x;

public:
  Lit() = default;
  constexpr Lit(const Lit &) = default;
  constexpr Lit(Var var, bool sign = false) : x(var + var + (int)sign) {}

  inline bool operator==(Lit p) const { return x == p.x; }
  inline bool operator!=(Lit p) const { return x != p.x; }

  inline bool operator<(Lit p) const {
    return x < p.x; // '<' makes p, ~p adjacent in the ordering.
  }

  inline Lit operator~() const {
    Lit p;
    p.x = x ^ 1;
    return p;
  }

  inline Lit operator^(bool b) const {
    Lit p;
    p.x = x ^ (unsigned int)b;
    return p;
  }

  inline bool sign() const { return x & 1; }
  inline int var() const { return x >> 1; }
  inline int toInt() const { return x; }
};

constexpr Lit lit_Undef(var_Undef, false);
constexpr Lit lit_Error(var_Undef, true);

//=================================================================================================
// Lifted booleans:
//
// NOTE: this implementation is optimized for the case when comparisons between values are mostly
//       between one variable and one constant. Some care had to be taken to make sure that gcc
//       does enough constant propagation to produce sensible code, and this appears to be somewhat
//       fragile unfortunately.

#define l_True  (lbool((uint8_t)0)) // gcc does not do constant propagation if these are real constants.
#define l_False (lbool((uint8_t)1))
#define l_Undef (lbool((uint8_t)2))

class lbool {
    uint8_t value;

public:
    explicit lbool(uint8_t v) : value(v) { }

    lbool()       : value(0) { }
    explicit lbool(bool x) : value(!x) { }

    bool  operator == (lbool b) const { return !!(((b.value&2) & (value&2)) | (!(b.value&2)&(value == b.value))); }
    bool  operator != (lbool b) const { return !(*this == b); }
    lbool operator ^  (bool  b) const { return lbool((uint8_t)(value^(uint8_t)b)); }

    bool  isTrue()  const { return *this == l_True;  }
    bool  isFalse() const { return *this == l_False; }

    lbool operator && (lbool b) const {
        uint8_t sel = (this->value << 1) | (b.value << 3);
        uint8_t v   = (0xF7F755F4 >> sel) & 3;
        return lbool(v); }

    lbool operator || (lbool b) const {
        uint8_t sel = (this->value << 1) | (b.value << 3);
        uint8_t v   = (0xFCFCF400 >> sel) & 3;
        return lbool(v); }

    friend int   toInt  (lbool l);
    friend lbool toLbool(int   v);
};
inline int   toInt  (lbool l) { return l.value; }
inline lbool toLbool(int   v) { return lbool((uint8_t)v);  }

//=================================================================================================
// Clause -- a simple class for representing a clause:

class Clause;
typedef Clause *CRef;
const CRef CRef_Undef = NULL;

class Clause {
private:
  std::vector<Lit> literals;
  float activity;

public:
  Clause(const vec<Lit> &ps, bool learnt) : literals(ps.size()) {
    for (int i = 0; i < ps.size(); i++)
      literals[i] = ps[i];

    activity = learnt ? 0.0f : -1.0f;
  }

  Clause(std::vector<Lit> &lits, bool learnt)
      : literals(std::move(lits)), activity(learnt ? 0.0f : -1.0f) {}

  int size() const { return literals.size(); }
  bool learnt() const { return activity >= 0.0f; }

  Lit &operator[](int i) { return literals[i]; }
  Lit operator[](int i) const { return literals[i]; }

  float get_activity() {
    assert(activity >= 0.0f);
    return activity;
  }

  void set_activity(float a) {
    assert(activity >= 0.0f && a >= 0.0f);
    activity = a;
  }
};

//=================================================================================================
// Watcher (moved here from Solver.h)

struct Watcher {
    CRef cref;
    Lit  blocker;
    Watcher(CRef cr, Lit p) : cref(cr), blocker(p) {}
    bool operator==(const Watcher& w) const { return cref == w.cref; }
    bool operator!=(const Watcher& w) const { return cref != w.cref; }
};

//=================================================================================================
}

#endif
