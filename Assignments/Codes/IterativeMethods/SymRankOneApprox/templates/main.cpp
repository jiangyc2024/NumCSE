#include "symrank1.hpp"

int main() {

  MatrixXd M(3, 3);
  double eps = 1e-4;
  M << 1 + eps, 2, 1, 2, 1 + eps, 1, 1, 1, 1 + eps;
  /*
   * run symRankOneBestApproxSym
   */
  VectorXd z = symRankOneBestApproxSym(M);

  std::cout << "The output of symRankOneBestApproxSym is:\n" << z << "\n";

  /*
   * test of symRankOneBestApproxSym
   */

  Vector3d ztest = { 1.21315, 1.21315,0.888086};

  if (ztest.size() == z.size()) {
    if ((ztest - z).norm() < eps)
      std::cout << "Test of symRankOneBestApproxSym passed!\n\n";
    else
      std::cout << "Test of symRankOneBestApproxSym failed: wrong output.\n\n";
  } else
    std::cout << "Test of symRankOneBestApproxSym failed: wrong size.\n\n";

  VectorXd v(3);
  v << 2, 0, 1;
  VectorXd b(9);
  b << 1, 2, 3, 2, 3, 4, 3, 4, 5;
  /*
   * run computeKronProdVecMult
   */

  VectorXd w = computeKronProdVecMult(v, b);

  std::cout << "The output of computeKronProdVecMult is:\n" << w << "\n";

  /*
   * test of computeKronProdVecMult
   */

  Vector3d wtest = {10, 16, 22};

  if (wtest.size() == w.size()) {
    if (wtest == w)
      std::cout << "Test of computeKronProdVecMult passed!\n\n";
    else
      std::cout << "Test of computeKronProdVecMult failed: wrong output.\n\n";
  } else
    std::cout << "Test of computeKronProdVecMult failed: wrong size.\n\n";

  /*
   * run symmRankOneApprox
   */
  MatrixXd MM(3, 3);
  MM << 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0;
  VectorXd zz = symmRankOneApprox(MM);

  std::cout << "The output of symmRankOneApprox is:\n" << zz << "\n";

  /*
   * test of symmRankOneApprox
   */

  Vector3d zztest = { 1.16657,0.441537,0.859084};

  if (zztest.size() == zz.size()) {
    if ((zztest - zz).norm() < eps)
      std::cout << "Test of symmRankOneApprox passed!\n\n";
    else
      std::cout << "Test of symmRankOneApprox failed: wrong output.\n\n";
  } else
    std::cout << "Test of symmRankOneApprox failed: wrong size.\n\n";
}
