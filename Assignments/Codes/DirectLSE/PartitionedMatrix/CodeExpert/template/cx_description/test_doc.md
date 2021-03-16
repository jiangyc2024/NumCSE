

# Tests of blockdecomp.hpp

> There are two tests.

***
> 1. We create a "random" 10 by 10 matrix of the given structure ($n=9$):
```
    int n = 9;
    std::srand(9);
    Eigen::MatrixXd R = Eigen::MatrixXd::Random(n,n).triangularView<Eigen::Upper>();
    Eigen::VectorXd v = Eigen::VectorXd::Random(n);
    Eigen::VectorXd u = Eigen::VectorXd::Random(n);
    Eigen::VectorXd bb = Eigen::VectorXd::Random(n+1);
```
Then we use solvelse() to solve the LSE and display the (transposed) solution.
The test is passed if the correct solution is printed:
```
-15.5807 -41.5599 -40.5345 -25.0211  4.18201  -8.0062 -1.49806  9.27245 -8.95513 -2.55923
```

***
> 2. The function
```
bool testSolveLSE(const MatrixXd & R, const VectorXd & v, const VectorXd & u, const VectorXd & b, VectorXd & x)
```
should return `true` if the solution from `solvelse()`, $y$, is close enough to the `Eigen::partialPivLU` solution, $x$. Moreover, it should save the `Eigen::partialPivLU` solution in the vector `x` which is passed as the last argument.
The second test is passed if `testSolveLSE()` returns `true` *and* $\| x - y \|$ is close to zero (tolerance $10^{-6}$).
