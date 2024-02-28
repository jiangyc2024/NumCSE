///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016-2020 SAM, D-MATH
/// Author(s): Ralf Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace QuadFormCOO {

using COOMat = std::vector<Eigen::Triplet<double>>;

Eigen::MatrixXd cooToDense(const COOMat &A) {
  Eigen::Index m = 0;
  Eigen::Index n = 0;
  for (const Eigen::Triplet<double> &t : A) {
    m = (m > t.row()) ? m : t.row() + 1;
    n = (n > t.row()) ? n : t.col() + 1;
  }
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(m, n);
  for (const Eigen::Triplet<double> &t : A) {
    M(t.row(), t.col()) += t.value();
  }
  return M;
}

/* SAM_LISTING_BEGIN_1 */
double evalQuadFormCOO(const COOMat &A, const Eigen::VectorXd &x) {
  double s = 0.0;
  for (const Eigen::Triplet<double> &t : A) {
    assert((t.row() < x.size()) && (t.col() < x.size()));
    s += (1.0 - x[t.row()]) * x[t.col()] * t.value();
  }
  return s;
}
/* SAM_LISTING_END_1 */

} // namespace QuadFormCOO

int main(int /*argc*/, char ** /*argv*/) {

  const QuadFormCOO::COOMat A{
      Eigen::Triplet<double>(0, 0, 1.0), Eigen::Triplet<double>(1, 1, 1.0),
      Eigen::Triplet<double>(2, 2, 1.0), Eigen::Triplet<double>(0, 1, 2.0),
      Eigen::Triplet<double>(1, 0, 2.0), Eigen::Triplet<double>(1, 2, 3.0),
      Eigen::Triplet<double>(2, 1, 3.0)};
  const Eigen::VectorXd x = (Eigen::VectorXd(3) << 1.0, 2.0, 3.0).finished();

  const Eigen::MatrixXd M{QuadFormCOO::cooToDense(A)};
  std::cout << "A = " << M << std::endl;
  std::cout << "x^T*A*x = " << QuadFormCOO::evalQuadFormCOO(A, x) << " <-> "
            << (Eigen::Vector3d::Ones() - x).transpose() * M * x << std::endl;
  return 0;
}
