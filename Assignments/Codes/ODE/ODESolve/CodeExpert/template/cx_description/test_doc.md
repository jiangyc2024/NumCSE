

# Tests of odesolve.hpp

***
> `Vector psitilde(DiscEvlOp &&Psi, unsigned int p, double h, const Vector &y0)`: is tested for the explicit Euler ($p = 1$) and equation $\dot{Y} = Y, Y(0) = 1$, with $h = 0.1$. One step returns 
```
1.105
```
> The same equation and `Psi` are used to test `std::vector<Vector> odeintequi(DiscEvlOp &&Psi, double T, const Vector &y0, unsigned int N)`. With `T = 1` and `N = 8` the result is
```
1
1.125
1.26562
1.42383
1.60181
1.80203
2.02729
2.2807
```

> The code `double testcvpExtrapolatedEuler(void)` generates the table
```
N       Error
4       0.0650468
8       0.0209801
16      0.00595508
32      0.00158261
64      0.000407501
128     0.000103355
256     2.60232e-05
512     6.52882e-06
1024    1.63508e-06
2048    4.09129e-07
4096    1.02327e-07
```
Only the second column is checked in the test. There is a separate test for the correct fitted rate
= 1.94797, obtained by polyfit.

***

> `std::pair <std::vector<double>, std::vector<Vector>> odeintssctrl(DiscEvlOp&& Psi, double T, const Vector &y0, double h0, unsigned int p, double reltol, double abstol, double hmin)`: the result is
```
1
1.01005
1.02444
1.03903
1.05382
1.06883
1.08405
1.09949
```
with the call `odeintssctrl(Psi, T, y0, 0.01, 1, 10e-5, 10e-5, 10e-5).second` (same equation as above and Euler method). There is no test for the plot.