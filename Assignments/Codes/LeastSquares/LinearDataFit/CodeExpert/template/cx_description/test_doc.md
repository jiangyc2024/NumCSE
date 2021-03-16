# Tests of linfit.hpp

> There are 5 functions in `linfit.hpp` and `main.cpp` contains the following tests:

***

> `MatrixXd make_A(const VectorXd &b)`: the expected result of this test is

```
10      100     2.4596   6.04965
5       25      2.22554  4.95303
3.33333 11.1111 2.01375  4.0552
2.5     6.25    1.82212  3.32012
2       4       1.64872  2.71828
1.66667 2.77778 1.49182  2.22554
1.42857 2.04082 1.34986  1.82212
1.25    1.562   1.2214   1.49182
1.11111 1.23457 1.10517  1.2214
1       1       1        1
```

containing the evaluations of the "span" functions for the linear fitting at the sample abscissae.

***
For the functions
> `VectorXd data_fit_normal(const MatrixXd &A, const VectorXd &b)`

> `VectorXd data_fit_qr(const MatrixXd &A, const VectorXd &b)` 

the students decide which one to test entering "1" or "2" respectively to the console. In both cases the correct result is the vector $\gamma$ that equals
```
4.05912 
0.614035 
-2.53146 
0.70576
```

***
The plots are also done with the method according to the  number "1" or "2" that the student decided. The results are very similar. There is no test for the following functions

> `void fit_plot( const VectorXd &gamma, const VectorXd &tl, VectorXd &yl)`
                
> `void error_plot(  const VectorXd &gamma, const MatrixXd &A, const VectorXd &f, VectorXd &err )`

In these functions the student should compute the values to plot for the fitting and for the error. The actual plot of these values is already done in the `main.cpp` using `matplotlibcpp`.