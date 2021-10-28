///
/// Minimal runner for (9-14)
///

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

#include "symrank1.hpp"

int main() {
  // Subtask (9-14.a)
  //  equ. (9.14.1) for a symmetric matrix M
  Eigen::MatrixXd M(3, 3);
  constexpr double eps = 1e-4;
  M << 1 + eps, 2, 1, 2, 1 + eps, 1, 1, 1, 1 + eps;
  /*
   * run symRankOneBestApproxSym
   */
  Eigen::VectorXd z = symRankOneBestApproxSym(M);

  std::cout << "The output of symRankOneBestApproxSym is:\n" << z << "\n";
  std::cout << "\n";

  // Subtask (9-14.f)
  //  efficient function to evaluate (9.14.14)
  Eigen::VectorXd v(3);
  v << 2, 0, 1;
  Eigen::VectorXd b(9);
  b << 1, 2, 3, 2, 3, 4, 3, 4, 5;
  /*
   * run computeKronProdVecMult
   */

  Eigen::VectorXd w = computeKronProdVecMult(v, b);

  std::cout << "The output of computeKronProdVecMult is:\n" << w << "\n";
  std::cout << "\n";

  // Subtask (9-14.g)
  //  Gauss-Newton iteration to solve equ. (9.14.4)
  Eigen::MatrixXd MM(3, 3);
  MM << 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0;
  /*
   * run symmRankOneApprox
   */
  Eigen::VectorXd zz = symmRankOneApprox(MM);

  std::cout << "The output of symmRankOneApprox is:\n" << zz << "\n";
}
