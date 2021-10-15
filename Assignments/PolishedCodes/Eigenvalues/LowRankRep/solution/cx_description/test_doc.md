# Tests of lowrankrep.hpp

***
> `std::pair<MatrixXd,MatrixXd> factorize_X_AB(const MatrixXd& X,const unsigned int k)`: We test this function on the matrix $X$ defined by
```
5 0
2 1
7 4
```
with $rank(X)=2$. Correct outputs are $A=$
```
4.62726 -1.89431
2.22977 0.167728
7.99362  1.04977 
```
and
$B=$
```
0.925453 -0.378863
0.378863  0.925453
```

Note that other combinations are possible (e.g. $AD$ and $BD^{-1}$ are solutions for any $D$ diagonal) therefore we only test the equality $X=AB^T$ up to the tolerance 1e-6.

***

> `std::tuple<MatrixXd, MatrixXd, MatrixXd> svd_AB(const MatrixXd& A, const MatrixXd& B)`: The input for the test are $A$
```
2 1 
2 3 
6 1
```
and
$B=$
```
4 4
5 0
7 4
3 2
```
The output, using correctly `Eigen::JacobiSVD` is $U=$
```
-0.321736 0.0871625
-0.458572  0.856829
-0.828371 -0.508179
```
$S=$
```
78.1139       0
      0 7.08661
```
and $V=$
```
-0.463767  0.557885
-0.418033 -0.819215
-0.714587 0.0663564
 -0.31549    0.1151
```
The three matrices are tested separately.

***

> `std::pair<MatrixXd,MatrixXd> rank_k_approx(const MatrixXd & Ax,const MatrixXd & Ay, const MatrixXd & Bx,const MatrixXd & By)`: Here the input for the test are $A_X, A_Y,B_X,B_Y$ given by
```
 0  9   8 -2   2  1    4  4
 2  6   3  4   2 -3   -5  0
 3  5   5  8   6  7    3  2
-2  3   7  5   8  4    5  9
 4  8   6  3   5  1    0  5
 9  0   2  1
```
respectively. Correct outputs are $A_Z=$
```
-111.401  57.2614
  -129.9 -3.24889
-187.161 -21.9279
-105.832  17.5224
 -178.86  6.35771
-117.068 -41.3817
```
and
$B_Z=$
```
-0.315078 0.0510129
  0.24722 -0.765217
-0.530351  0.253255
 -0.70135 -0.296917
-0.257768 -0.509454
```
Again multiple solutions are correct so we test only $A_Z(B_Z^T) = Z$ where $Z=$
```
 38.021  -71.3579   73.5833   61.1292 -0.456398
40.7628  -29.6277   68.0697   92.0698   35.1392
57.8517  -29.4904   93.7078   137.776   59.4155
34.2391  -39.5721   60.5656   69.0223   18.3532
56.6791  -49.0827   96.4687   123.556   42.8655
34.7744   2.72452   51.6068   94.3923   51.2584
```
