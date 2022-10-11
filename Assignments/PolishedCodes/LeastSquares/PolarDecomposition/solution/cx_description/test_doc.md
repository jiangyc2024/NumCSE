# Tests of polardecomposition.hpp

> `void initialize` is tested using $m=5$, $n=4$ on a random matrix. The decomposition is initialized, the original matrix is restored and the symmetry and correctness are checked. Finally, we test if $Q$ has orthogonal columns.

> `PolarDecomposition()` is tested using $m=11$, $n=7$ and $k=3$ on random matrices. The same checks as with the constructor are performed.