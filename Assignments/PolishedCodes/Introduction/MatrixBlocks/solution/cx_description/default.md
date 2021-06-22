
# Tutorial: part 2/3
## Matrix blocks


> Eigen matrices have a range of methods to access their submatrices (blocks), instead of individual entries at a time. Suppose A is an Eigen::MatrixXd object, then a few useful methods are:
* `A.row(i)`: the i-th row of A,
* `A.col(j)`: the j-th column of A,
* `A.topRows(p)`: the first p rows of A,
* `A.leftCols(q)`: the first q columns of A,
* `A.bottomRightCorner(p,q)`: the p by q matrix formed by the intersection of the last p rows and last q columns of A, and more generally
* `A.block(i,j,p,q)`: the p by q submatrix of A, whose top-left corner is at the index (i,j).


> The purpose of this exercise is to learn to declare and initialize Eigen matrices, and to get familiar with some handy typedefs and methods.

> Open "MatrixBlocks.hpp" and fill in the missing code in between the delimiters `// START` and `// END` according to the instructions preceded by `// TODO:`.

> Read each function of the program carefully. The code will not compile until after the first two tasks are finished, which is to include headers and set the namespace.