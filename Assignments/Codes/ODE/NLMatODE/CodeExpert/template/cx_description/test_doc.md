

# Tests of NLmatode.hpp

***
> `MatrixXd matode(const MatrixXd &Y0, double T)`: We test this function on the matrix $Y_0$ defined by
```
1  1  0
0  3  2
1  5  2
```
and $T=1$. The result $Y(1)$ should be 
```
  1.11717   1.06478 -0.317651
  0.85854   5.24544   2.58701
 0.121829   2.52023   1.09839
```

***
> `bool checkinvariant(const MatrixXd &M, double T)`: is tested for 2 different matrices: $Y_0$ above and $M_0 =$
```
1  1  0
1  3  1
0  1  1
```
$M_0$ is symmetric, so the equation is stationary and the invariant is preserved. 
The invariant is not preserved for $Y_0$.  

***
> `double cvgDiscreteGradientMethod(void)`: runs the method with $Y_0 =$
```
0 1 0 0 0
0 0 1 0 0
0 0 0 1 0
0 0 0 0 1
1 0 0 0 0
```
and the resulting rate is 1.98098. It also prints the table
```
M       Error
10      0.00260457
20      0.000645272
40      0.000160726
80      4.00747e-05
160     9.9619e-06
320     2.44364e-06
640     5.84909e-07
1280    1.98872e-07
```
which is not tested.