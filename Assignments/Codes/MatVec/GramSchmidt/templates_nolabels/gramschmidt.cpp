#include <Eigen/Dense>
#include <iostream>

using Eigen;

/* @brief Performs Gram-Schidt orthonormalization
 * Given a matrix $\mathbf{A}$ of linearly independent columns,
 * returns the result of a Gram-Schmidt orthonormalization.
 * Ustable GS algorithm: output is prone to cancellation issues.
 * @param[in] $\mathbf{A}$ Matrix of linearly independent columns
 * @return Matrix with ONB of $span(a_1, \cdots, a_n)$ as columns
 */
MatrixXd gram_schmidt(const MatrixXd & A) {
  // We create a matrix Q with the same size as A
  Matrix Q(A);


  return Q;
}

int main(void) {
  // Orthonormality test
  unsigned int n = 9;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
  Eigen::MatrixXd Q = gramschmidt( A );

  // Norm of matrices?
  double err = (Q*Q.transpose() - Identity(n,n)).norm();

  // Output should be identity matrix
  std::cout << "Error: "
            << err
            << std::endl;

  double eps = std::numeric_limits<double>::denorm_min();
  exit(err < eps);
}
