

# Tests of shape_ident.hpp

> There are three functions in shape_ident.hpp and main.cpp has one test for each function.
***
> `shape_ident_matrix( MatrixXd X )`: This function returns the matrix B.
The test prints `shape_ident_matrix(Xstop)`, and the correct result is:
```
 1 -3  0  0
 0  0  1 -3
 3 -1  0  0
 0  0  3 -1
 3  1  0  0
 0  0  3  1
 1  3  0  0
 0  0  1  3
-1  3  0  0
 0  0 -1  3
-3  1  0  0
 0  0 -3  1
-3 -1  0  0
 0  0 -3 -1
-1 -3  0  0
 0  0 -1 -3
```
***
> `solve_lsq( MatrixXd X, MatrixXd P, MatrixXd A )`: This function calculates the least squares solution, saves it in A and returns the norm of the residual.
The test prints `solve_lsq( Xstop, P1, A )` and `A`, and the correct result is:
```
7.47635
0.0202649 -0.0374318
0.00595625  0.0872423

```
***
> `solve_lsq( MatrixXd Xstop, MatrixXd Xpriority, MatrixXd P, MatrixXd A )`: This function decides whether points P are more similar to Xstop or Xpriority.
The test prints `solve_lsq( Xstop, Xpriority, P1, A )`, `solve_lsq( Xstop, Xpriority, P2, A )` and `solve_lsq( Xstop, Xpriority, P3, A )`, and the correct result is:
```
Priority
Stop
Stop
```
