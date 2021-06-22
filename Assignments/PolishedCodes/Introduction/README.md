# Problems

## MatrixClass (1-1)

- the solution of 1-1.a is technically not exactly correct. Including <Eigen/Core> would be sufficient to just use matrices. Dense includes the following: Core, LU, Cholesky, QR, SVD, Geometry, Eigenvalues
- changed interface of constantTriangular: use unsigned ints for indexing