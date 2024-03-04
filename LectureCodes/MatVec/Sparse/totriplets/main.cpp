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
  trp.emplace_back(0, 0, 3.0);
  trp.emplace_back(0, 2, 4.0);
  trp.emplace_back(1, 2, 1.0);
  trp.emplace_back(2, 1, 2.0);
  trp.emplace_back(2, 3, 5.0);
  trp.emplace_back(2, 3, 1.0);
  trp.emplace_back(3, 0, 4.0);

  const int rows = 4;
  const int cols = 4;
  Eigen::SparseMatrix<double> A(rows, cols);
  A.setFromTriplets(trp.begin(), trp.end());
  A.makeCompressed();
  std::cout << "Sparse matrix A = " << std::endl << A << std::endl;
  // Extract triplets
  const std::vector<Eigen::Triplet<double>> extr_trip{convertToTriplets(A)};
  // List triplets
  for (const Eigen::Triplet<double> &triplet : extr_trip) {
    std::cout << "triplet (" << triplet.row() << ", " << triplet.col() << ", "
              << triplet.value() << ")" << std::endl;
  }
  return 0;
}
