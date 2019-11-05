

# Tests of adaptivepolyintp.hpp

***
> ` VectorXd adaptivepolyintp( Function&& f, double a, double b, double tol, int N )` : 
We call this function with data $a = 0, b = 1, tol = 10^{-4}, N = 100$ on the function $f_1(t) = sin(e^{2t})$. The vector of interpolation nodes should be
```
      0.5
 0.777778
        1
        0
 0.191919
 0.929293
0.0707071
 0.636364
 0.333333
 0.979798
 0.030303
 0.858586
 0.262626
```

***

> `VectorXd adaptivepolyintp( Function&& f, double a, double b, double tol, int N, std::vector<double> *errortab = nullptr )`: 
It is is called with same data as before, and returns also the vector `errortab`. The test is passed if the following values are returned.
```
n                \eps_n
-------------------------
1               1.41046
2               3.0219
3               12.6848
4               1.13185
5               0.453926
6               0.423054
7               0.148235
8               0.0463045
9               0.0173202
10              0.0147524
11              0.00483056
12              0.00128332
```

***

> `void plotInterpolationError(void)`: this a function produces plots for $n$ vs errors for the input data $tol= 10^{-6}, N=1000$ and for the functions $f_1$ above as well as $f_2(t) = \frac{\sqrt{t}}{1 + 16 t^2}$. There is no test for the plot.
