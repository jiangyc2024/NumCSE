

# Tests of distfitting.hpp

> There are four functions in distfitting.hpp and main.cpp has one test for each function plus an additional test for the runtime of `solveNormalEquations()`.
***
> `initA( n )`: This function returns the matrix A.
The test accepts any result that is a row-permutation of:
```
-1., 1., 0.,
-1., 0., 1.,
-1., 0., 0.,
 0., -1., 1.,
 0., -1., 0.,
 0., 0., -1.;
```
***
> `solveExtendedNormalEquations( D )`: This function solves the least-squares problem described by the parameters in D using the Extended Normal Equations.
The test prints the resulting vector x:
```
 2
-1
-2
```
***
> `solveNormalEquations( D )`: This function solves the least-squares problem described by the parameters in D using the Sherman-Morrison-Woodbury formula.
The test prints the resulting vector x:
```
 2
-1
-2
```
It also tests if the runtime is approximately in O(n^2) by calling `solveNormalEquations( D )` with matrices D of size 2^L for $L = 3,...,7$. A timeout might occur if your implementation is to slow.
***
> `testNormalEquations( D )`: This functions tests whether `solveExtendedNormalEquations()` and `solveNormalEquations()` return approximately the same result.
The test simply checks whether the function returns true.
