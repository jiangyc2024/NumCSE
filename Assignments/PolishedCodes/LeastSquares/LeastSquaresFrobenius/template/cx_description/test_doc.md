# Tests of matrixfit.hpp

***
> `Eigen::MatrixXd min_frob(const VectorXd & z, const VectorXd & g)`: We test this function on the vectors defined by
```
z = [1,2,3,4,5]^T
g = [5,4,3,2,1]^T
```
and the result is tested with the matrix $M^*=$
```
0.0909091  0.181818  0.272727  0.363636  0.454545
0.0727273  0.145455  0.218182  0.290909  0.363636
0.0545455  0.109091  0.163636  0.218182  0.272727
0.0363636 0.0727273  0.109091  0.145455  0.181818
0.0181818 0.0363636 0.0545455 0.0727273 0.0909091
```
***

> `min_frob` is additionally tested in `tests.cpp` against the solution.