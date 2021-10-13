#include <Eigen/Dense>
#include <iostream>

#include "distfitting.hpp"
#include "timer.h"

int main() {
  /*
   *	run initA
   */
  constexpr unsigned int n = 4;
  Eigen::MatrixXd Adense = Eigen::MatrixXd(initA(n));

  std::cout << "A(4) is:\n" << Adense << std::endl;
  /*
   *	test of initA
   */
  Eigen::MatrixXd Atest(6, 3), Atest2(6, 3);
  // clang-format off
  Atest << -1., 1., 0., 
           -1., 0., 1., 
           -1., 0., 0., 
            0., -1., 1.,
            0., -1., 0.,
            0., 0., -1.;
  // clang-format on

  if (Adense.size() == Atest.size()) {
    // Test for each row of Adense if it exists in Atest
    const unsigned int n = Adense.rows();
    bool all_rows_contained = true;  // True if Adense is a permutation of Atest

    for (unsigned int i = 0; i < n; ++i) {  // Loop over all rows of Adense
      Eigen::VectorXd row = Adense.row(i).transpose();
      bool contained = false;  // True if row i is contained in Atest
      for (unsigned int j = 0; j < n; ++j) {  // Loop over all rows of Atest
        Eigen::VectorXd row_test = Atest.row(j).transpose();
        if ((row - row_test).norm() == 0) contained = true;
      }
      if (contained == false) all_rows_contained = false;
    }

    if (all_rows_contained)
      std::cout << "Test for initA passed!\n\n";
    else
      std::cout << "Test for initA failed: wrong output.\n\n";
  } else
    std::cout << "Test for initA failed: wrong size.\n\n";

  /*
   *	run solveExtendedNormalEquations
   */
  Eigen::MatrixXd D(4, 4);
  D << 0.0, -3.0, -4.0, -2.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0,
      0.0, 0.0;

  Eigen::VectorXd x = solveExtendedNormalEquations(D);
  std::cout << "Solution extended normal equations:\n" << x << std::endl;
  /*
   *	test of solveExtendedNormalEquations
   */
  Eigen::VectorXd xtest(3);
  xtest << 2, -1, -2;

  if (x.size() == xtest.size()) {
    if ((x - xtest).norm() < 1e-10)
      std::cout << "Test for solveExtendedNormalEquations passed!\n\n";
    else
      std::cout
          << "Test for solveExtendedNormalEquations failed: wrong output.\n\n";
  } else
    std::cout
        << "Test for solveExtendedNormalEquations failed: wrong size.\n\n";

  /*
   *	run solveNormalEquations
   */
  Eigen::VectorXd y = solveNormalEquations(D);
  std::cout << "Solution normal equations:\n" << y << std::endl;
  /*
   *	test of solveNormalEquations
   */
  if (y.size() == xtest.size()) {
    if ((y - xtest).norm() < 1e-10)
      std::cout << "Test for solveNormalEquations passed!\n\n";
    else
      std::cout << "Test for solveNormalEquations failed: wrong output.\n\n";
  } else
    std::cout << "Test for solveNormalEquations failed: wrong size.\n\n";
  /*
   * test of runtime of solveNormaEquations
   */

  constexpr unsigned int steps = 5;  // testing until dim = 2^(steps+2)
  Eigen::VectorXd time(steps);

  for (unsigned int s = 0; s < steps; ++s) {
    const unsigned int dim = std::pow(2, s + 3);
    Eigen::MatrixXd rand = Eigen::MatrixXd::Zero(dim, dim);
    rand.triangularView<Eigen::Upper>() = Eigen::MatrixXd::Random(dim, dim);
    rand.diagonal() = Eigen::VectorXd::Zero(dim);

    Timer timer = Timer();
    for (unsigned int i = 0; i < 10; ++i) {
      timer.start();
      Eigen::VectorXd y_test = solveNormalEquations(rand);
      timer.stop();
    }
    time(s) = timer.mean();
  }

  Eigen::VectorXd time_log = time.array().log().matrix();
  Eigen::VectorXd coeff = time_log.tail(steps - 1) - time_log.head(steps - 1);
  coeff = coeff / std::log(2);
  std::cout << "Your runtime is O(n^" << coeff.maxCoeff() << ")" << std::endl;
  if (coeff.maxCoeff() < 2.5)
    std::cout << "Test for runtime of solveNormalEquations passed!\n\n";
  else
    std::cout
        << "Test for runtime of solveNormalEquations failed: too slow.\n\n";

  const bool same = testNormalEquations(D);
  if (same)
    std::cout << "Test for testNormalEquations passed!\n\n";
  else
    std::cout << "Test for testNormalEquations failed.\n\n";
}