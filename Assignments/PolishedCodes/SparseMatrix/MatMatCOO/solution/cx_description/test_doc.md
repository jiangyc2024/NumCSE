# Tests of matmatCOO.hpp

> main.cpp has one test for each function of gridfun.hpp.

***
> `void eval(MatrixXd & X, std::function<double(index_t,index_t)> f))`: The function should create the matrix $X$ defined by
```
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```
of size $(n,m) = (16,16)$.
***
> `SparseMatrix<double> build_matrix(const Matrix3d & S, const shape_t & size)`: We test this function with size $(m,n) = (3,3)$, as in 3-16.b. The result is:

```
-4  1  0  1  0  0  0  0  0 
1 -4  1  0  1  0  0  0  0 
0  1  -4  0  0  1  0  0  0 
1  0  0 -4  1  0  1  0  0 
0  1  0  1 -4  1  0  1  0 
0  0  1  0  1 -4  0  0  1 
0  0  0  1  0  0 -4  1  0  
0  0  0  0  1  0  1 -4  1 
0  0  0  0  0  1  0  1 -4
```
***

> `void mult(const SparseMatrix<double> & A,
const MatrixXd & X, MatrixXd & Y)`: multiplication of $A*vec(X)$ where $X$ is the matrix above and the matrix $A$ is given by
`build_matrix` with input size $(16,16)$. The exact solution is $Y$ given by 
```
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  1  1  1  1  1  1  1  0  0  0  0  0
 0  0  0  1 -2 -1 -1 -1 -1 -1 -2  1  0  0  0  0
 0  0  0  1 -1  0  0  0  0  0 -1  1  0  0  0  0
 0  0  0  1 -1  0  0  0  0  0 -1  1  0  0  0  0
 0  0  0  1 -1  0  0  0  0  0 -1  1  0  0  0  0
 0  0  0  1 -1  0  0  0  0  0 -1  1  0  0  0  0
 0  0  0  1 -1  0  0  0  0  0 -1  1  0  0  0  0
 0  0  0  1 -2 -1 -1 -1 -1 -1 -2  1  0  0  0  0
 0  0  0  0  1  1  1  1  1  1  1  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 ```
 
 ***
 
 > `void solve(const SparseMatrix<double> & A,
const MatrixXd & Y, MatrixXd & X)`: stores in $X$ the solution of the system $A*vec(X) = vec(Y)$ where $Y$ is the matrix above. The exact solution is again $X$: here we test the Frobenius norm of the solution, which exact value equals $7$

> In tests.cpp every function to be implemented is tested against the solution.