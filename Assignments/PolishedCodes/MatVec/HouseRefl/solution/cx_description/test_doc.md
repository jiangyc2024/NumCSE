# Tests of houserefl.hpp

> There is one function in `houserefl.hpp`:
* `void houserefl(const Eigen::VectorXd &v, Eigen::MatrixXd &Z)`

> This function is tested in `main.cpp` for a random vector `v` of dimension
> $n=6$. The test utilizes the fact that a for a orthogonal matrix we have
> $Z^T Z =  Z Z^T = Id$.

> The test is passed if the error is less than $10^{-10}$.
