

# Tests of multAmin.hpp

> There are four functions in multAmin.hpp:
* `multAminSlow(const VectorXd & x, VectorXd & y)` (not written by the students),
* `multAmin(const VectorXd & x, VectorXd & y)`,
* `multAmin_runtime(void)` (not tested since it only prints runtimes), and
* `multABunitv()`.

> `multAmin(const VectorXd & x, VectorXd & y)` is tested with
```
unsigned int M = 10;
VectorXd xa = VectorXd::Random(M);
VectorXd ys, yf;
multAmin(xa, yf);
multAminSlow(xa, ys);
```
The error `(ys - yf).norm()` is printed, and the test is passed if the result is less that $10^{-6}$.

> `multABunitv()` is just tested by calling `multABunitv()`, the test is only passed if the student writes the printing command themselves. The program should print
```
C = 
1 0 0 0 0 0 0 0 0 0
0 1 0 0 0 0 0 0 0 0
0 0 1 0 0 0 0 0 0 0
0 0 0 1 0 0 0 0 0 0
0 0 0 0 1 0 0 0 0 0
0 0 0 0 0 1 0 0 0 0
0 0 0 0 0 0 1 0 0 0
0 0 0 0 0 0 0 1 0 0
0 0 0 0 0 0 0 0 1 0
0 0 0 0 0 0 0 0 0 1
```
as well as return that matrix.
