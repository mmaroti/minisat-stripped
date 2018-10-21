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
#include <unordered_map>

#include "minisat/mtl/Alg.h"
#include "minisat/mtl/Vec.h"
#include "minisat/mtl/Map.h"

namespace Minisat {

//=================================================================================================
// Variables, literals, lifted booleans, clauses:


// NOTE! Variables are just integers. No abstraction here. They should be chosen from 0..N,
// so that they can be used as array indices.

typedef int Var;
#define var_Undef (-1)


struct Lit {
    int     x;

    // Use this as a constructor:
    friend Lit mkLit(Var var, bool sign);

    bool operator == (Lit p) const { return x == p.x; }
    bool operator != (Lit p) const { return x != p.x; }
    bool operator <  (Lit p) const { return x < p.x;  } // '<' makes p, ~p adjacent in the ordering.
};


inline  Lit  mkLit     (Var var, bool sign = false) { Lit p; p.x = var + var + (int)sign; return p; }
inline  Lit  operator ~(Lit p)              { Lit q; q.x = p.x ^ 1; return q; }
inline  Lit  operator ^(Lit p, bool b)      { Lit q; q.x = p.x ^ (unsigned int)b; return q; }
inline  bool sign      (Lit p)              { return p.x & 1; }
inline  int  var       (Lit p)              { return p.x >> 1; }

// Mapping Literals to and from compact integers suitable for array indexing:
inline  int  toInt     (Var v)              { return v; }
inline  int  toInt     (Lit p)              { return p.x; }
inline  Lit  toLit     (int i)              { Lit p; p.x = i; return p; }

//const Lit lit_Undef = mkLit(var_Undef, false);  // }- Useful special constants.
//const Lit lit_Error = mkLit(var_Undef, true );  // }

const Lit lit_Undef = { -2 };  // }- Useful special constants.
const Lit lit_Error = { -1 };  // }


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

class Clause {
    struct {
        unsigned mark      : 2;
        unsigned learnt    : 1;
        unsigned size      : 27; }                            header;
    union data_union { Lit lit; float act; } data[0];

    friend class ClauseAllocator;

    // NOTE: This constructor cannot be used directly (doesn't allocate enough memory).
    template<class V>
    Clause(const V& ps, bool learnt) {
        header.mark      = 0;
        header.learnt    = learnt;
        header.size      = ps.size();

        for (int i = 0; i < ps.size(); i++)
            data[i].lit = ps[i];

        if (header.learnt)
            data[header.size].act = 0;
    }

public:
    int          size        ()      const   { return header.size; }
    void         shrink      (int i)         { assert(i <= size()); if (header.learnt) data[header.size-i] = data[header.size]; header.size -= i; }
    void         pop         ()              { shrink(1); }
    bool         learnt      ()      const   { return header.learnt; }
    uint32_t     mark        ()      const   { return header.mark; }
    void         mark        (uint32_t m)    { header.mark = m; }
    const Lit&   last        ()      const   { return data[header.size-1].lit; }

    // NOTE: somewhat unsafe to change the clause in-place! Must manually call 'calcAbstraction' afterwards for
    //       subsumption operations to behave correctly.
    Lit&         operator [] (int i)         { return data[i].lit; }
    Lit          operator [] (int i) const   { return data[i].lit; }
    operator const Lit* (void) const         { return (Lit*)data; }

    float&       activity    ()              { assert(header.learnt); return data[header.size].act; }
};


//=================================================================================================
// ClauseAllocator -- a simple class for allocating memory for clauses:


const CRef CRef_Undef = NULL;
class ClauseAllocator 
{
    static int clauseWord32Size(int size, bool learnt){
        return (sizeof(Clause) + (sizeof(Lit) * (size + (int)learnt))) / sizeof(uint32_t); }
 public:
    template<class Lits>
    CRef alloc(const Lits& ps, bool learnt = false)
    {
        assert(sizeof(Lit)      == sizeof(uint32_t));
        assert(sizeof(float)    == sizeof(uint32_t));

        Clause *cid = (Clause*)new uint32_t[clauseWord32Size(ps.size(), learnt)];
        new (cid) Clause(ps, learnt);

        return cid;
    }

    void free(CRef r) { 
        delete[]((uint32_t*)r);
    }

    // Deref, Load Effective Address (LEA), Inverse of LEA (AEL):
    Clause*       lea       (CRef r)       { return (Clause*)r; }
    const Clause* lea       (CRef r) const { return (Clause*)r; }
};


//=================================================================================================
// OccLists -- a class for maintaining occurence lists with lazy deletion:

template<class Idx, class Vec, class Deleted>
class OccLists
{
    vec<Vec>  occs;
    vec<char> dirty;
    vec<Idx>  dirties;
    Deleted   deleted;

 public:
    OccLists(const Deleted& d) : deleted(d) {}

    void  init      (const Idx& idx){ occs.growTo(toInt(idx)+1); dirty.growTo(toInt(idx)+1, 0); }
    // Vec&  operator[](const Idx& idx){ return occs[toInt(idx)]; }
    Vec&  operator[](const Idx& idx){ return occs[toInt(idx)]; }
    Vec&  lookup    (const Idx& idx){ if (dirty[toInt(idx)]) clean(idx); return occs[toInt(idx)]; }

    void  cleanAll  ();
    void  clean     (const Idx& idx);
    void  smudge    (const Idx& idx){
        if (dirty[toInt(idx)] == 0){
            dirty[toInt(idx)] = 1;
            dirties.push(idx);
        }
    }

    void  clear(bool free = true){
        occs   .clear(free);
        dirty  .clear(free);
        dirties.clear(free);
    }
};


template<class Idx, class Vec, class Deleted>
void OccLists<Idx,Vec,Deleted>::cleanAll()
{
    for (auto const& d : dirties) {
        if (dirty[toInt(d)]) {
            clean(d);
        }
    }
    dirties.clear();
}


template<class Idx, class Vec, class Deleted>
void OccLists<Idx,Vec,Deleted>::clean(const Idx& idx)
{
    Vec& vec = occs[toInt(idx)];
    auto j = std::remove_if(vec.begin(), vec.end(), [&] (typename Vec::value_type const& elem) {
            return deleted(elem);
        }
    );
    vec.truncate(j);
    dirty[toInt(idx)] = 0;
}


//=================================================================================================
// CMap -- a class for mapping clauses to values:


template<class T>
class CMap
{
    struct CRefHash {
        uint32_t operator()(CRef cr) const { return (uint32_t)cr; } };

    std::unordered_map<CRef, T, CRefHash> map;

 public:
    // Size-operations:
    void     clear       ()                           { map.clear(); }
    int      size        ()                const      { return map.elems(); }


    // Insert/Remove/Test mapping:
    void     insert      (CRef cr, const T& t){ map.insert(cr, t); }
    void     growTo      (CRef cr, const T& t){ map.insert(cr, t); } // NOTE: for compatibility
    void     remove      (CRef cr)            { map.remove(cr); }
    bool     has         (CRef cr, T& t)      { return map.peek(cr, t); }

    // Vector interface (the clause 'c' must already exist):
    const T& operator [] (CRef cr) const      { return map[cr]; }
    T&       operator [] (CRef cr)            { return map[cr]; }
};

//=================================================================================================
}

#endif
