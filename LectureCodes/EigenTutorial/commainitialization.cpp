// **********************************************************************
// Eigen demo code for course "Numerical Methods for CSE"
// Problem of comma initialization in Eigen
// Contributed by Yann Billeter
// See commainitialization_notes.pdf for explanations
// **********************************************************************

#include <Eigen/Dense>
#include <iostream>

int main() {
  std::cout << "Demonstration of different results obtained by comma "
               "initializatoion and"
            << std::endl
            << "direct assignment in Eigen" << std::endl;
  std::cout << "(Contributed by Yann Billeter, October 2021)" << std::endl;
  // system dimension
  const unsigned int n = 9;
  // test matrix and vector
  Eigen::MatrixXd R(n, n);
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      R(i, j) = 1.0 / (1.0 + i + j);
    }
  }
  Eigen::VectorXd b = Eigen::VectorXd::LinSpaced(n, 1.0, n);
  Eigen::VectorXd xa = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd xb = Eigen::VectorXd::Zero(n);

  //----------------------------------------------------------------
  xa << R.triangularView<Eigen::Upper>().solve(b);
  xb = R.triangularView<Eigen::Upper>().solve(b);
  //----------------------------------------------------------------

  std::cout << "LSE solution A: " << xa.transpose() << std::endl;
  std::cout << "LSE solution B: " << xb.transpose() << std::endl;

  std::cout << "Error compared between A and B: " << (xa - xb).norm()
            << std::endl;
}
