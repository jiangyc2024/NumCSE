#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iomanip>
#include <iostream>
#include <vector>

#include "tridiagleastsquares.hpp"

int main() {
  std::cout << "Enter problem dimension (n>=4): ";
  std::size_t n;
  std::cin >> n;
  // Define tridiagonal matrix.
  constexpr double alpha = -4.0;
  constexpr double beta = 1.0;
  Eigen::SparseMatrix<double> T(n, n);
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(3 * n - 2);
  for (std::size_t i = 0; i < n - 1; ++i) {
    triplets.push_back(Eigen::Triplet<double>(i, i, alpha));
    triplets.push_back(Eigen::Triplet<double>(i + 1, i, beta));
    triplets.push_back(Eigen::Triplet<double>(i, i + 1, beta));
  }
  triplets.push_back(Eigen::Triplet<double>(n - 1, n - 1, alpha));
  T.setFromTriplets(triplets.begin(), triplets.end());
  T.makeCompressed();

  // Synthetic data
  Eigen::VectorXd z = Eigen::VectorXd::LinSpaced(n, -1.0, 1.0);
  std::srand(39);
  Eigen::VectorXd c = T * z + 0.5 * Eigen::VectorXd::Random(n);
  Eigen::VectorXd x = lsqEst(z, c);

  std::cout << "Least squares parameters: \n"
            << "(alpha, beta) = " << std::setprecision(6) << x.transpose()
            << std::endl;
}