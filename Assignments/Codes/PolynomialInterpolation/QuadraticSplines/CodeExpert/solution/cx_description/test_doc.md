# Tests of quadsplines.hpp

> Note that you can use functions from subproblems that you did not implement correctly. The tests will call these from the solution. However, every function will be tested.

> `std::pair<VectorXd, VectorXd> increments(const VectorXd &t)`: We test this function on the values $t = (.1,.2,.5)$  (computations are pointless for equispaced nodes) and the result should be
```
0.5
0.1
0.1
0.3
0.5
0.1
0.6
0.2
0.4
0.8
0.6
```

***

> `VectorXd compute_c( const VectorXd &t, const VectorXd &y)`: using $y = (1,2,1,2 )$ and the times above, the correct result is $c=$
```
0.355736
 1.15991
0.209733
 1.40556
```
 
 ***
 
 > `VectorXd compute_d ( const VectorXd &c, const VectorXd &t )`: computes the values of the spline at the points $(0,t_1,...,t_{n-1},1)$, that is
 ```
 1.06141
1.51564
1.84473
1.31634
1.06141
 ```
***

> `VectorXd quadspline(const VectorXd &t, const VectorXd &y, const VectorXd &x)`: we evaluate the spline computed with $t$ and $y$ above at the points $x = (0,.4,.9)$. Note that at 0 we recover $d_0$ (consistency check). The result is
```
1.06141
0.976439
1.63152
```

***

> `void plotquadspline( const std::string &filename )`: there is no test for the plot. 

***

> `std::vector<double> qsp_error(unsigned int q)`: should print the table
```
n       Error values: 
------------------------ 
2       0.606776
4       0.434326
8       0.0402719
16      0.00253016
32      0.000261929
64      3.13621e-05
128     3.87463e-06
256     4.8303e-07
```

