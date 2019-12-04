

# Tests of matode.hpp

***

We test the three methods on $h = 0.1$, the matrix $A =$
```
0 1 0
1 0 1
1 1 0
```
and initial value $Y_0 = Id$.
> The results after one step are: `MatrixXd eeulstep(const MatrixXd & A, const MatrixXd & Y0, double h)`
```
  1   0.1     0
0.1     1   0.1
0.1   0.1     1
```
`MatrixXd ieulstep(const MatrixXd & A, const MatrixXd & Y0, double h)`
```
1.01124   0.102145  0.0102145
0.11236   1.02145   0.102145
0.11236   0.11236   1.01124
```
and `MatrixXd impstep(const MatrixXd & A, const MatrixXd & Y0, double h)`
```
 1.00528   0.100515  0.00502576
 0.105541  1.0103    0.100515
 0.105541  0.105541  1.00528
```

***

> `std::tuple<double,double,double> checkOrthogonality(void)`: the student defines the input $A =$
```
0  1  1
-1 0  1
-1 -1 0
```
and $Y_0$ as the matrix $Q$ of the QR decomposition of $M =$
```
8  1  6
3  5  7
9  9  2
```
The 20th step is tested: it should return the values
```
0.00850951     0.00845861    9.21619e-15
```
for the three methods.