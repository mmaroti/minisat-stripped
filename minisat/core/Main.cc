
#include <iostream>
#include <vector>

#include "minisat/core/Solver.h"

using namespace Minisat;

int main() {
  Solver solver;

  // binary relation
  const int size = 8;
  Lit table[size][size];
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      table[i][j] = mkLit(solver.newVar());

  // reflexive
  for (int i = 0; i < size; i++)
    solver.addClause(table[i][i]);

  // symmetric
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      solver.addClause(~table[i][j], table[j][i]);

  // transitive
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      for (int k = 0; k < size; k++)
        solver.addClause(~table[i][j], ~table[j][k], table[i][k]);

  // find all solutions
  int solutions = 0;
  while (solver.solve()) {
    solutions += 1;
    std::vector<Lit> clause(size * size);
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++) {
        bool b = solver.modelValue(table[i][j]).isTrue();
        clause[i * size + j] = table[i][j] ^ b;
      }
    solver.addClause(clause);
  }

  std::cout << "solutions: " << solutions;
  if (solutions != 4140)
    std::cout << " is incorrect";
  std::cout << std::endl;

  return 0;
}
