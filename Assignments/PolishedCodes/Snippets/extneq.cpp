/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: October 2022
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Extneq {
/* SAM_LISTING_BEGIN_1 */
using Triplet = Eigen::Triplet<double>;
using TripletMatrix = std::vector<Triplet>;
TripletMatrix extNeqSysMatCOO(unsigned int m, unsigned int /*n*/,
                              const TripletMatrix &A_coo) {
  TripletMatrix A_ext_coo{};
  // Set up the $m\times m$ identity matrix block
  for (int j = 0; j < (int)m; ++j) {
    A_ext_coo.push_back(Triplet{j, j, -1.0});
  }
  for (const auto &t : A_coo) {
    // Right upper block $\cob{\VA}$
    A_ext_coo.push_back(Triplet{t.row(), t.col() + (int)m, t.value()});
    // Left lower block $\cob{\VA^\top}$
    A_ext_coo.push_back(Triplet{t.col() + (int)m, t.row(), t.value()});
  }
  return A_ext_coo;
}
/* SAM_LISTING_END_1 */
}  // namespace Extneq

int main(int /*argc*/, char ** /*argv*/) {
  Extneq::TripletMatrix A_coo{{0, 0, 2.0}, {0, 1, 3.0}, {1, 0, 4.0}};
  auto A_ext_coo = Extneq::extNeqSysMatCOO(2, 2, A_coo);
  Eigen::SparseMatrix<double> A_ext(4, 4);
  A_ext.setFromTriplets(A_ext_coo.begin(), A_ext_coo.end());
  Eigen::MatrixXd A = A_ext;
  std::cout << "A = " << std::endl << A << std::endl;
  return 0;
}
