

# Tests of kron.hpp

> There are four functions in kron.hpp:
* `kron(const MatrixXd &A, const MatrixXd &B, MatrixXd &C)`
* `kron_mult(const MatrixXd &A, const MatrixXd &B, const VectorXd &x, VectorXd &y)`
* `kron_reshape(const MatrixXd &A, const MatrixXd &B, const VectorXd &x, VectorXd &y)`
* `kron_runtime()`

> The last one `kron_runtime()` prints the runtimes of the other three functions. There is no test for this.

> The other three functions are tested on identical inputs, and should all return the same output.
The inputs are defined as
```
  MatrixXd A(2, 2);
  A << 1, 2, 3, 4;
  MatrixXd B(2, 2);
  B << 5, 6, 7, 8;
  srand(5);
  VectorXd x = Eigen::VectorXd::Random(4);
```
and the functions should output
```
y =
-7.919
-10.05
-23.53
-30.51
```