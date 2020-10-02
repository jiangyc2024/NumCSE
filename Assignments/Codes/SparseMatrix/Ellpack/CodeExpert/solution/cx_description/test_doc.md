

# Tests of ellpack.hpp

> There are two functions and a construntor in ellpack.hpp: main.cpp has one test for each.

***
> `EllpackMat(const Triplets & triplets, index_t m, index_t n);`: In the output of main, We print the result of this constructor on the matrix E defined by the following triplets
```
(1, 2, 4)
(0, 0, 5)
(1, 2, 6)
(2, 5, 7)
(0, 4, 8)
(0, 0, 1)
(1, 3, 9)
(2, 2, 10)
(1, 3, 2)
(2, 1, 11)
(1, 0, 12)
```

***

> `mvmult()`: it should compute $Ex$, where
$x = [4,5,6,7,8,9]^{T}$. The result is compared with the exact solution  $[88,185,178]^{T}$.

> `mtvmult()`: it should compute $E^Tx$ where $x=[1,2,3]^T$. The result is compared with the exact solution $[30,33,50,22,8,21]^{T}$.

>Addition to the output from main, We test the correctness of constructor, maxcols, `mvmult()` and `mtvmult()` using random matrix of size $70 * 100$ with 0.1 repeating ratio. 
>For multiplication, we also use random generated vectors of length 100 and 70 respectively
>and check the correctness of their multiplication with the random matrix and its tranpose respectively.