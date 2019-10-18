

# Tests of matrixfit.hpp

***
> `MatrixXd min_frob(const VectorXd & z, const VectorXd & g)`: We test this function on the vectors defined by
```
z = [1,2,3,4,5]^T
g = [5,4,3,2,1]^T
```
and the result is tested with the matrix $M^*=$
```
0.0909091  0.181818  0.272727  0.363636  0.454545
0.0727273  0.145455  0.218182  0.290909  0.363636
0.0545455  0.109091  0.163636  0.218182  0.272727
0.0363636 0.0727273  0.109091  0.145455  0.181818
0.0181818 0.0363636 0.0545455 0.0727273 0.0909091
```
***

> `bool testMformula(unsigned int n)`: This function is a test written by the student, so there is no explicit test in Code Expert. $n$ is the input size of the vectors $z$ and $g$. We run it with $n=100$. The students are asked to define the vectors
```
z = [1,2,...,100]^T
g = [1,-1,1,...,-1,1,-1]^T
```
The test is passed if the exact formula returns the same output as the function above, up to tolerance 1e-10.