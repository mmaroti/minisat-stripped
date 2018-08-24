# Notes about the inputs

The basic set of inputs were taken from the [SATLIB benchmark page]
(http://www.cs.ubc.ca/~hoos/SATLIB/benchm.html). For simplicity, they are
split into two main folders, `SAT` and `UNSAT`, depending on whether the
problem is satisfiable or unsatisfiable.

The `easy.txt` file contains a list of inputs that finish in reasonable
time and thus can be run per commit in a CI setting.

### DIMACS Benchmark Instances (original source: [DIMACS Benchmark set for SAT](ftp://dimacs.rutgers.edu/pub/challenge/satisfiability/benchmarks/cnf/))
* AIM: Artificially generated Random-3-SAT - 48 instances satisfiable, 24 unsatisfiable
* LRAN: Large Random-3-SAT instances - 3 instances, all satisfiable
* JNH: Random SAT instances with variable length clauses - 16 instances satisfiable, 34 instances unsatisfiable
* DUBOIS: Randomly generated SAT instances - 13 instances, all unsatisfiable
* GCP: Large SAT-encoded Graph Colouring problems - 4 instances, all satisfiable
* PARITY: Instances for problem in learning the parity function - 20 instances, all satisfiable
* II: Instances from a problem in inductive inference - 41 instances, all satisfiable
* HANOI: SAT-encoding of Towers of Hanoi - 2 instances, all satisfiable
* BF: Circuit fault analysis: bridge fault - 4 instances, all unsatisfiable
* SSA: Circuit fault analysis: single-stuck-at fault - 4 instances satisfiable, 4 instances unsatisfiable
* PHOLE: Pigeon hole problem - 5 instances, all unsatisfiable
* PRET: Encoded 2-colouring forced to be unsatisfiable - 8 instances, all unsatisfiable


