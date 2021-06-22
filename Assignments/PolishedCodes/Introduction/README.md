# Problems

## MatrixClass (1-1)

- the solution of 1-1.a is technically not exactly correct. Including <Eigen/Core> would be sufficient to just use matrices. Dense includes the following: Core, LU, Cholesky, QR, SVD, Geometry, Eigenvalues
- changed interface of constantTriangular: use unsigned ints for indexing

## MatrixBlocks (1-2)

- removed 
````
// We can use the Eigen namespace to improve the readability of our code.
// This allows us to skip the "Eigen::" in "Eigen::MatrixXd" for example.
// TO DO: Add the following line below: using namespace Eigen;
// START
using namespace Eigen;
// END
````
because omitting the namespace does not comply with Google coding guidelines
- small cosmetic changes