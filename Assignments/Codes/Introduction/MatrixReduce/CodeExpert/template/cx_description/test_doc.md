

# Tests of MatrixReduction.hpp

> There are four functions in MatrixReduction.hpp and main.cpp has one test for each function.

***

> `average( MatrixXd A )`: This function should calculate the average of the entries in A.
The test prints `average(Matrix3d::Identity())`, and the correct result is:
```
  0.333333
```
***

> `percent_zero( MatrixXd A )`: This function calculates the percentage of entries in A that are exactly equal to zero.
The test prints `percent_zero(Matrix3d::Identity())`, and the correct result is:
```
  66.6667
```
***
> `has_zero_column( MatrixXd A )`: This function checks if any one of the columns in A is equal to zero.
The test prints 
`has_zero_column(B)` followed by `has_zero_column(B.transpose())` where `B=MatrixXd::Identity(4,5)`, and the correct result is:
```
  1 0
```
***
> `columns_sum_to_zero( MatrixXd A )`: This functions returns a copy of A, except that the diagonal entries have been changed in a way such that the columns of the new matrix sum to zero.
The test prints `columns_sum_to_zero(C)` where `C=Matrix3d::Random()+Matrix3d::Constant(1)` (the random seed is fixed), and the correct result is:
```
-0.721194  0.160061  0.561133
0.0929355 -0.989846  0.896911
  1.98551  0.159289   -2.1448
```
