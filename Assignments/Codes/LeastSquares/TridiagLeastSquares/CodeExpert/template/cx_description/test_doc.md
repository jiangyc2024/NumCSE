

# Tests of tridiagleastsquares.hpp


> To test `linReg()` we define `T` as the $999\times 999$ tridiagonal matrix with diagonal entries equal to $\alpha=-4$, and off-diagonal entries equal to $\beta=1$.
Then we define $z$ and $c$ according to
```
  int n=99
  VectorXd z = VectorXd::LinSpaced(n,-1.0,1.0);
  std::srand(39);
  VectorXd c = T*z + 0.5*VectorXd::Random(n);
```
The test is passed if `lsqEst(z,c)` has the correct output
$x = (-4.1839,\ 1.09028)^T$.

