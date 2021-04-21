# Tests of legendre.hpp

> `legvals` is tested against the solution with with a linspaced vector of size $N = 20$ and matrices of size $20 \times 10$.

> There is no test for `Pnx` as it is an optional helper function.

> `gaussPts` is tested against the solution with $n = 8$ and default tolerances. It will pass even though `gaussPts` gives technically wrong results. However, it is tested whether the implementation matches the solution's implementation.

> `gaussPts_regulaFalsi` is tested against the solution with $n = 8$ and default tolerances.
