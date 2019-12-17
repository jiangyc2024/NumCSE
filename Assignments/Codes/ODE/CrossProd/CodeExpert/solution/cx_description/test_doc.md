

# Tests of cross.hpp

***
> `std::vector<VectorXd> solve_imp_mid(Function &&f, Jacobian &&Jf, double T, const VectorXd &y0, unsigned int N)`: this is tested on the ODE $\dot{y} = [y_1y_2, y_2y_3, y_3 - y_1]^T$ with initial data $y(0) = [0.1,0.2, 0.4]^T$ and parameters: T = 1, N = 1. The result is
```
0.135782
0.407031
0.964218
```
The same test is repeated for the function 
`std::vector<VectorXd> solve_lin_mid(Function &&f, Jacobian &&Jf, double T, const VectorXd &y0, unsigned int N)`. With the results
```
0.131724
0.371034
0.968276
```