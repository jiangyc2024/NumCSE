

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
```
The output, using correctly `Eigen::JacobiSVD` is $U=$
```
-0.319939 -0.0935435
-0.441442  -0.865779
-0.838313   0.491605
```
$S=$
```
48.781       0
      0 6.95785
```
and $V=$
```
-0.740879 -0.671638
-0.671638  0.740879
```
The three matrices are tested separately.

***

> `std::pair<MatrixXd,MatrixXd> rank_k_approx(const MatrixXd & Ax,const MatrixXd & Ay, const MatrixXd & Bx,const MatrixXd & By)`: Here the input for the test are $A_X, A_Y,B_X,B_Y$ given by
```
1 0       8 -2      2  1      4  4
9 2       3  4      2  3      5  0
6 3       5  8    
```
respectively. Correct outputs are $A_Z=$
```
-46.9539  15.3405
-61.8361 -1.13729
-80.8728 -8.03698
```
and
$B_Z=$
```
-0.764387 -0.644758
-0.644758  0.764387
```
Again multiple solutions are correct so we test only $A_ZB_Z = Z$ where $Z=$
```
26 42
48 39
67 46
```
