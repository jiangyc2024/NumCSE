# Tests of matrixlowrank.hpp

> For the tests of `MatrixLowRank::operator*()`, `MatrixLowRank::operator=*()`
> and `MatrixLowRank::addTo()` we define the following test data:
```
  double tol = 1e-4;
  MatrixXd A(4,2),B(3,2),C(3,2),D(4,2),E(3,1),F(4,1),G(3,4);
  A << 1,2,3,4,5,6,7,8;
  B << 9,0,1,2,3,4;
  C << 5,6,7,8,9,0;
  D << 179,94,437,250,695,406,953,562;
  E << 9.+1e-9,1+1e-9,3+1e-9;
  F << -1,-3,-5,-7;  
  G << 0,0,0,0,4,8,12,16,8,16,24,32;
  MatrixLowRank L(B,A);
  MatrixLowRank M(A,B);
  MatrixLowRank N(E,F);
```

> To test `MatrixLowRank::operator*()` we compute `MC = M * C` which should be
> equal to `D`.

> To test `MatrixLowRank::operator=*()` we compute `M *= C` which should be
> equal to `D`.

> To test `MatrixLowRank::addTo()` we compute `N.addTo(L)` as well as the
> SVD of `G`. We then compare the rank of the `SVD` of G and the rank of `N`.


> There are 5 functions in `linfit.hpp` and `main.cpp` contains the following tests:
