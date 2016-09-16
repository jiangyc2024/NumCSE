///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Xiaolin Guo, Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

# include <iostream>
# include "./denseSolve.hpp"
# include "./sparseSolve.hpp"

void test(const MatrixXd& A, const VectorXd& b) {
  VectorXd xLu, xQr, xSvd, xSp, xSpd;
  lu_solve(A, b, xLu);
  qr_solve(A, b, xQr);
  svd_solve(A, b, xSvd);
  sparse_solve(A.sparseView(), b, xSp);

  // build SPD matrix (strongly diagonally dominant & hermitian)
  auto upper = A.template triangularView<Eigen::Upper>();
  MatrixXd A_spd = upper.toDenseMatrix() + upper.toDenseMatrix().transpose();
  // A_spd now is hermitian, lets make it strongly diagonally dominant
  for (int i = 0; i < A.rows(); ++i) {
    A_spd(i,i) = 2*(A_spd.row(i).cwiseAbs().sum());
  }

  
  sparseSpd_solve(A_spd.sparseView(), b, xSpd);

  std::cout << "-- Error w/ size n = " << b.size() << " ------------\n"
            << "LU:     " << (A*xLu - b).norm() << "\n"
            << "QR:     " << (A*xQr - b).norm() << "\n"
            << "SVD:    " << (A*xSvd - b).norm() << "\n"
            << "SPD:    " << (A_spd*xSpd - b).norm() << "\n"
            << "Sparse: " << (A*xSp - b).norm() << "\n";
}

int main() {
  for (int n = 5; n < 50; n *= 2) {
    test(Eigen::MatrixXd::Random(n,n), Eigen::VectorXd::Random(n));
  }

  return 0;
}
