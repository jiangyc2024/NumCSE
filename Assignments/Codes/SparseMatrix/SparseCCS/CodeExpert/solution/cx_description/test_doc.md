

# Tests of sparseCCS.hpp

> 

***
> `void CCS(const MatrixXd & A, VectorXd & val, VectorXd & row_ind, VectorXd & col_ptr)`: We test this function on the matrix A defined by

```
4, -1,  0, -1,  0,  0,
-1,  4, -1,  0, -1, 0,
0, -1,  4,  0,  0, -1,
-1,  0,  0,  4, -1, 0,
0, -1,  0, -1,  4, -1,
0,  0, -1,  0, -1,  4;

```
> The tests are passed if the l^2 error of the vectors val, row_ind and col_ptr are < 10^{-6}. The reference values are computed with the built-in

`SparseMatrix<double> As = A.sparseView();
As.makeCompressed();`