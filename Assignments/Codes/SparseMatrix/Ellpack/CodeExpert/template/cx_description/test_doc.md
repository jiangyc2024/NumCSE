

# Tests of ellpack.hpp

> There are two functions and a construntor in ellpack.hpp: main.cpp has one test for each.

***
> `EllpackMat(const Triplets & triplets, index_t m, index_t n);`: We test this constructor on the matrix E defined by the following triplets
```
(1, 2, 4)
(0, 0, 5)
(1, 2, 6)
(2, 5, 7)
(0, 4, 8)
(1, 3, 9)
(2, 2, 10)
(2, 1, 11)
(1, 0, 12)
```

***

> `mvmult()`: it should compute $Ex$, where
$x = [4,5,6,7,8,9]^{T}$. The result is compared with the exact solution  $[84,171,178]^{T}$.

> `mtvmult()`: it should compute $E^Tx$ where $x=[1,2,3]^T$. The result is compared with the exact solution $[29,33,50,18,8,21]^{T}$.