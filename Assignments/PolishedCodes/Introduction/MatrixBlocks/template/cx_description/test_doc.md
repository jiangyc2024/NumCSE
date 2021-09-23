

# Tests of MatrixBlocks.hpp

> There are three functions in MatrixBlocks.hpp and there is one test for each function.
***
> `zero_row_col( MatrixXd A, int p, int q )`: This function returns a copy of A with row p and column q set to zero.
The test prints `zero_row_col(Matrix3d::Constant(-1),0,1)`, and the correct result is:
```
 0  0  0
-1  0 -1
-1  0 -1
```
***
> `swap_left_right_blocks( MatrixXd A, int p )`: This function swaps splits the matrix A into two blocks, the first p columns and the rest, and returns a matrix in which these blocks have been swapped.
The test prints `swap_left_right_blocks(MatrixXd::Identity(4,3),2)`, and the correct result is:
```
0 1 0
0 0 1
1 0 0
0 0 0
```
***
> `tridiagonal( int n, double a, double b, double c)`: This function creates an n by n tridiagonal matrix with constant diagonals with the values a, b, and c.
The test prints `tridiagonal(4,-1,2,-1)`, and the correct result is:
```
 2 -1  0  0
-1  2 -1  0
 0 -1  2 -1
 0  0 -1  2
```
