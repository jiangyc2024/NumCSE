///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s):
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Sparse>
#include <iostream>

Eigen::SparseMatrix<double, Eigen::RowMajor> reserve_demo(unsigned int n) {
  const Eigen::Index rows = n;
  const auto cols = 2 * static_cast<Eigen::Index>(n);
  const unsigned int max_no_nnz_per_row = 3;
  Eigen::SparseMatrix<double, Eigen::RowMajor> mat(rows, cols);
  mat.reserve(Eigen::RowVectorXi::Constant(rows, max_no_nnz_per_row));
  // do many (incremental) initializations
  for (Eigen::Index i = 0; i < rows; ++i) {
    mat.insert(i, i) = -1.0;  // only for matrix entries not yet set!
    mat.insert(i, i + 1) = 1.0;
    mat.coeffRef(i, 2 * i) -= 1.0;  // access entry possibly not set yet
  }
  mat.makeCompressed();  // squeeze out zeros
  return mat;
}

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "C++ demo program for course Numerical Methods for CSE"
            << std::endl;
  const Eigen::SparseMatrix<double, Eigen::RowMajor> mat{reserve_demo(21)};
  const Eigen::MatrixXd M = mat;
  std::cout << "M = " << M.rows() << " x " << M.cols()
            << ", outerSize = " << mat.outerSize()
            << ", nnz = " << mat.nonZeros() << ", comp = " << mat.isCompressed()
            << " : \n " << M << std::endl;
  return 0;
}
