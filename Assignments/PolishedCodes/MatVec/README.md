IMPORTANT: Problems marked with ! didn't had a CodeExpert folder i.e. were
changed more than the others. Filenames might have changed which might break
latex.

# Issues in Latex
- Most problems don't have a link to the template/github like e.g. problem
  1-4 but then maybe the issue is with problem 1-4 having such a link?
- Problem 1.7.b refers to code 1.7.5 but should be 1.7.2 I think? The normal
  problem collection does it correct.
- Problem 1-5.b has the word "gather" in the description. Not sure that makes
  sense.

# Problems
- Changed all 2-* to 1-* and 20 to 21 (as in 2020 to 2021) in all problems.

## ArrowMatrix (1-1)
- Basically deleted everything except CodeExpert folder
- Changed `arrow_matrix_2_times_x(d,a,x,yi);` to 
  `arrow_matrix_2_times_x(d,a,x,yi);` to times in `main.cpp`.

## GramSchmidt (1-2)
- Removed `zzzzz_test_runner.cpp`.
- Made usage of `Eigen` namespace explicit i.e. used `Eigen::`

## Kroenecker (1-3)
- No major changes. Just simple cleanup.

## FastMatMult (1-4)
- No major changes. Just copied the solution and template from the already
  existing CodeExpert folder and cleaned up the formatting a bit.

## HouseRefl (1-5) !
- Split the code from `houserefl.cpp` into `main.cpp` and `houserefl.hpp`.
- The latex labels are now inside `houserefl.hpp`.
- There is only one label and it starts at 1. Should it start at 0?
- There were no tests implemented. Old `main.cpp` only outputted results. I
  added a simple test to `main.cpp` which consists of checking the frobenius
  norm of `(Z.transpose*Z - Id).norm() << 10^-10`

## StructuredMatrixVector (1-7)
- Made usage of `Eigen` namespace explicit i.e. used `Eigen::`
- Changed return type of `multAmin_runtime()` from int to void
- Improved terminal output of tests in `main.cpp`
- General imrovements to comments and code structure

## Cancellation (1-8)
- No major changes. Just simple cleanup.

## ComplexRoot (1-11)
- No major changes. Just copied the solution and template from the already
  existing CodeExpert folder and cleaned up the formatting a bit.
