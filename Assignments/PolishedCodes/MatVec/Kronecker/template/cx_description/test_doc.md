# Tests of kron.hpp

> There are four functions in kron.hpp:
* `kron(const MatrixXd &A, const MatrixXd &B, MatrixXd &C)`
* `kron_mult(const MatrixXd &A, const MatrixXd &B, const VectorXd &x, VectorXd &y)`
* `kron_reshape(const MatrixXd &A, const MatrixXd &B, const VectorXd &x, VectorXd &y)`
* `kron_runtime()`

> The last one `kron_runtime()` prints the runtimes of the other three functions. There is no test for this.

> The other three functions are tested against the master solution.