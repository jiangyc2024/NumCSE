

# Tests of MatrixClass.hpp

> There are four functions in MatrixClass.hpp and main.cpp has one test for each function.
***
> `smallTriangular(double a, double b, double c)`: This function does not need to be edited, so the first test is passed automatically as long as the code compiles. The function takes as arguments three doubles and creates a 2 by 2 upper triangular matrix from those numbers.
The test prints `smallTriangular(1,2,3)`, and the correct result is:
```
  1 2
  0 3
```
***
> `constantTriangular(int n, double val)`: This function creates an n by n triangular MatrixXd with the entries in upper triangular part equal to `val`.
The test prints `constantTriangular(3,20)`, and the correct result is:
```
  20 20 20
   0 20 20
   0  0 20
```
***
> `casting()`: This function calculates the inner product of the integer vector `[2 1]` and the complex vector `[1-i, 5+i]` and returns the real part as a double.
The test prints `casting()`, and the correct result is:
```
  6
```
***
> `arithmetics(int m, int n)`: This function performs some arbitrary arithmetics on complex n by n matrices.
The test prints `arithmetics(4)`, and the correct result is:
```
(-41.4706,27.1176)
(-23.7692,42.8462)
(-6.47692,48.4154)
        (7.2,48.4)
```