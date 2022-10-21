# Problem 0-3: Eigen's reduction operations on matrices

> For the task description of this exercise, please also refer to [NCSEFL_Problems.pdf](
https://www.sam.math.ethz.ch/~grsam/NumMeth/HOMEWORK/NCSEFL_Problems.pdf). 

> Eigen matrices have a range of methods that reduce them to a single number. Suppose A is an Eigen::MatrixXd object, then a few useful methods are:
* `A.size()`: the number of entries in A,
* `A.sum()`: the sum of all entries in A,
* `A.prod()`: the product of all entries in A,
* `A.minCoeff()`: the minimum value amongst the entries of A,
* `A.maxCoeff()`: the maximum value amongst the entries of A, and
* `A.norm()`: The (Frobenius) norm of the matrix A.

> The [Eigen::Array class](https://eigen.tuxfamily.org/dox/group__TutorialArrayClass.html) is also featured in this exercise to perform coefficient-wise comparison.

> Open "MatrixReduce.hpp" and fill in the missing code in between the delimiters `// START` and `// END` according to the instructions preceded by `// TODO:`.

> Read each function of the program carefully. The code will not compile until after the first task is finished, which is to include headers.