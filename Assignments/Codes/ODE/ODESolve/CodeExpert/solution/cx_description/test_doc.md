

# Tests of gramschmidt.hpp

> There are two functions in gramschmidt.hpp and main.cpp has one test for each function.

***
> `gram_schmidt(const MatrixXd & A)`: The code generates the table
```
4       0.0650468
8       0.0209801
16      0.00595508
32      0.00158261
64      0.000407501
128     0.000103355
256     2.60232e-05
512     6.52882e-06
1024    1.63508e-06
2048    4.09129e-07
4096    1.02327e-07
```


***

> `orthogonality_test()`: This function is written by the student. It should create a random matrix A, calculate Q = gram_schmidt(A) and return err = norm(Q^TQ - Id).
The test is passed if err < 10^{-9}.