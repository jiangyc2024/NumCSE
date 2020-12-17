
#include "symrank1.hpp"

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace Eigen;

int main() {
  // Subtask (9-14.a)
  //  equ. (9.14.1) for a symmetric matrix M
  MatrixXd M(3, 3);
  double eps = 1e-4;
  M << 1 + eps, 2, 1, 2, 1 + eps, 1, 1, 1, 1 + eps;
  /*
   * run symRankOneBestApproxSym
   */
  VectorXd z = symRankOneBestApproxSym(M);

  std::cout << "The output of symRankOneBestApproxSym is:\n" << z << "\n";
  std::cout << "\n";

  // Subtask (9-14.f)
  //  efficient function to evaluate (9.14.14)
  VectorXd v(3);
  v << 2, 0, 1;
  VectorXd b(9);
  b << 1, 2, 3, 2, 3, 4, 3, 4, 5;
  /*
   * run computeKronProdVecMult
   */

  VectorXd w = computeKronProdVecMult(v, b);

  std::cout << "The output of computeKronProdVecMult is:\n" << w << "\n";
  std::cout << "\n";

  // Subtask (9-14.g)
  //  Gauss-Newton iteration to solve equ. (9.14.4)
  MatrixXd MM(3, 3);
  MM << 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0;
  /*
   * run symmRankOneApprox
   */
  VectorXd zz = symmRankOneApprox(MM);

  std::cout << "The output of symmRankOneApprox is:\n" << zz << "\n";
}
