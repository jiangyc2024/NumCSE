
# Tests of eval_deriv.hpp

> `evaldp()`: Test is passed if the correct value of $(p(x),p'(x))$ is printed for $p(x) = -x^4 - 4x^3 + 2x^2 + 9x + 4$,  $x = -2.2$.
The correct output is
```
  13.0464, -15.288
```

> `evaldp_naive()`: same as `evaldp()`.

> `polyTestTime()`: Test passed if "test passed: 1" is printed.

> `dipoleval()`: Test is passed if the derivative of interpolating polynomial of the data $(t_j,y_j),\ j=1,2,3,4$ defined by
```
    int n = 4;
    std::srand(199);
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(n,-1,1);
    Eigen::VectorXd y = Eigen::VectorXd::Random(n);
```
is evaluated at $x=(-2, -1,  0,  1,  2)$ and the correct result is printed.
The correct output is 
```
  -7.37389  -2.99816 0.0951584   1.90606   2.43454
```

> `dipoleval_alt()`: same as `dipoleval()`.

> `testDipolEval()`: Test passed if "test passed: 1" is printed.
