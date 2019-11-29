# Tests of quasilin.hpp

***
> `Eigen::VectorXd fixed_point_step(const Eigen::VectorXd &xk, const Eigen::VectorXd &b)`

and 

> `Eigen::VectorXd newton_step(const Eigen::VectorXd &x, const Eigen::VectorXd &b)`
We test these functions on the rhs defined by $b = [1,1,1,1,1]^T$ and current step $x = [\frac{1}{\sqrt{5}},\frac{1}{\sqrt{5}},\frac{1}{\sqrt{5}},\frac{1}{\sqrt{5}},\frac{1}{\sqrt{5}}]^T$. the output for the two methods is
```
0.211538
0.153846
0.173077
0.153846
0.211538
```
and
```
0.259273
0.188562
0.212132
0.188562
0.259273
```
respectively.

***
> `Eigen::VectorXd solveQuasiLinSystem(double rtol, double atol, Eigen::VectorXd &b)` and  `Eigen::VectorXd solveQLSystem_Newton(double rtol, double atol, Eigen::VectorXd &b)` are tested with the same input rhs $b$ as before. The tests expects the result
```
0.2422138484
0.1627846781
0.1951187913 
0.1627846781
0.2422138484
```