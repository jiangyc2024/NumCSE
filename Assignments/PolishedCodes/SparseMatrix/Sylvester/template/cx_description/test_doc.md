

# Tests of Sylvester.hpp

> The three functions are tested in `main.cpp` at your disposal.

> `solveDiagSylvesterEq()` is tested for
$
A = diag([1, 2, 3, 4, 5, 6, 7, 8]) \in \mathbb{R}^{8,8}
$
The test is passed if the output is `X = `
```
 0.5    0    0    0    0    0    0    0
   0  0.4    0    0    0    0    0    0
   0    0  0.3    0    0    0    0    0
   0    0    0 0.24    0    0    0    0
   0    0    0    0 0.19    0    0    0
   0    0    0    0    0 0.16    0    0
   0    0    0    0    0    0 0.14    0
   0    0    0    0    0    0    0 0.12
```

***
> `sparseKron()` is tested for the sparse s.p.d. matrix `A = `
```
    1.4   -0.11       0       0       0       0  -0.072       0
  -0.11       1       0  -0.027       0       0       0       0
      0       0       1   -0.07       0       0       0 -0.0082
      0  -0.027   -0.07     1.5       0   -0.12       0       0
      0       0       0       0       1   -0.12       0       0
      0       0       0   -0.12   -0.12     1.2       0 -0.0047
 -0.072       0       0       0       0       0     1.1   0.048
      0       0 -0.0082       0       0 -0.0047   0.048       1
```
The output is compared to `kron()` from a previous exercise. The test is passed if the difference between the output of `kron()` and `sparseKron()` is smaller than the tolerance $10^{-6}$.

***

> `solveSpecialSylvesterEq()` is tested for the same `A` as `sparseKron()`.
The test is passed if the output is `X = `
```
    0.47    0.007 -2.8e-05 -0.00031 -3.8e-06 -4.4e-05   0.0055  0.00054
   0.007      0.5 -0.00029    0.002 -6.1e-05 -0.00037  -0.0013  7.2e-05
-2.8e-05 -0.00029      0.5   0.0054 -0.00015 -0.00092  8.7e-05  1.5e-05
-0.00031    0.002   0.0054     0.46  -0.0016    0.012 -2.3e-05 -0.00015
-3.8e-06 -6.1e-05 -0.00015  -0.0016      0.5   0.0049  5.4e-06 -0.00011
-4.4e-05 -0.00037 -0.00092    0.012   0.0049     0.49  3.9e-05  0.00019
  0.0055  -0.0013  8.7e-05 -2.3e-05  5.4e-06  3.9e-05      0.5 -0.00098
 0.00054  7.2e-05  1.5e-05 -0.00015 -0.00011  0.00019 -0.00098      0.5
```

> Additionally, we provide tests against the solution of these functions.