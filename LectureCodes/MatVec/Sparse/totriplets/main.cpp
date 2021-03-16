///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2019 SAM, D-MATH
/// Author(s): R. Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include "totriplets.hpp"

int main(int /*argc*/, char ** /*argv*/) {
  // Initialize a (small) sparse matrix
  std::vector<Eigen::Triplet<double>> trp;
  trp.push_back(Eigen::Triplet<double>(0, 0, 3.0));
  trp.push_back(Eigen::Triplet<double>(0, 2, 4.0));
  trp.push_back(Eigen::Triplet<double>(1, 2, 1.0));
  trp.push_back(Eigen::Triplet<double>(2, 1, 2.0));
  trp.push_back(Eigen::Triplet<double>(2, 3, 5.0));
  trp.push_back(Eigen::Triplet<double>(2, 3, 1.0));
  trp.push_back(Eigen::Triplet<double>(3, 0, 4.0));

  int rows = 4, cols = 4;
  Eigen::SparseMatrix<double> A(rows, cols);
  A.setFromTriplets(trp.begin(), trp.end());
  A.makeCompressed();
  std::cout << "Sparse matrix A = " << std::endl << A << std::endl;
  // Extract triplets
  std::vector<Eigen::Triplet<double>> extr_trip{convertToTriplets(A)};
  // List triplets
  for (Eigen::Triplet<double> &triplet : extr_trip) {
    std::cout << "triplet (" << triplet.row() << ", " << triplet.col() << ", "
              << triplet.value() << ")" << std::endl;
  }
  return 0;
}
