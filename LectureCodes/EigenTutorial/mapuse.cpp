// **********************************************************************
// Eigen demo code for course "Numerical Methods for CSE"
// **********************************************************************

#include <Eigen/Dense>
#include <cassert>
#include <iostream>

// Extract even and odd components of a vector of even length
void evenoddsplit(const Eigen::VectorXd &v) {
  const Eigen::Index n = v.size();
  assert((n % 2 == 0) || "Vector length must be even");
  const Eigen::Index m = n / 2;
  // Even/odd splitting by merely reinterpreting data in memory
  const Eigen::VectorXd ve{
      Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned, Eigen::Stride<2, 2>>(
          v.data(), m, 1)};
  const Eigen::VectorXd vo{
      Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned, Eigen::Stride<2, 2>>(
          v.data() + 1, m, 1)};
  std::cout << "vector v = [" << v.transpose() << "]" << std::endl;
  std::cout << "even components = [" << ve.transpose() << "]" << std::endl;
  std::cout << "odd components = [" << vo.transpose() << "]" << std::endl;
}

template <typename MATRIX>
void vectorizedemo(const MATRIX &M) {
  // The type of the matrix entries
  using entry_t = typename MATRIX::Scalar;
  // Size of the matrix
  int n_rows = M.rows();
  int n_cols = M.cols();
  if (M.IsRowMajor) {
    // For a row major matrix the task is more difficult, we have to use a
    // suitable stride.
    const auto v = Eigen::Map<const Eigen::Matrix<entry_t, Eigen::Dynamic, 1>,
                              Eigen::Unaligned,
                              Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>{
        M.data(), n_rows * n_cols,
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(1, n_cols)};
    std::cout << "M r.m.: vec(M) = " << v.transpose() << std::endl;
  } else {
    // Little has to be done in the case of a row major matrix, because we can
    // just wrap the contiguous array of matrix entries into a column vector.
    const auto v = Eigen::Map<const Eigen::Matrix<entry_t, Eigen::Dynamic, 1>>{
        M.data(), n_rows * n_cols, 1};
    std::cout << "M c.m.: vec(M) = " << v.transpose() << std::endl;
  }
}

int main(int /*argc*/, char ** /*argv */) {
  std::cout << "Demo for use of the Map class in Eigen" << std::endl;
  {
    std::cout << "Demo I: even-odd splitting" << std::endl;
    std::vector<double> vd{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    Eigen::VectorXd v{Eigen::Map<Eigen::VectorXd>(vd.data(), vd.size(), 1)};
    evenoddsplit(v);
  }
  {
    std::cout << "Demo II: vectorization of a matrix" << std::endl;
    Eigen::MatrixXd M =
        (Eigen::MatrixXd(3, 4) << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
            .finished();
    std::cout << "Matrix = " << std::endl << M << std::endl;
    vectorizedemo(M);
  }
  {
    std::cout << "Demo II: vectorization of a matrix" << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> M =
        (Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(
             3, 4)
             << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12).finished();
    std::cout << "Matrix = " << std::endl << M << std::endl;
    vectorizedemo(M);
  }
  return 0;
}
