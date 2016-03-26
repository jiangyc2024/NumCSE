# include "../../DividedDifferences/Eigen/divdiff.hpp"

using Eigen::VectorXd;

/* Evaluation of polynomial in Newton basis (divided differences)
 * IN:  t = nodes (mutually different)
 *      y = values in t
 *      x = evaluation points (as Eigen::Vector)
 * OUT: p = values in x                                           */
void evaldivdiff(const VectorXd& t, const VectorXd& y, const VectorXd& x, VectorXd& p) {
  const unsigned n = y.size() - 1;

  // get coefficients of polynomial
  VectorXd coeffs;
  divdiff(t, y, coeffs);

  // evaluate
  VectorXd ones = VectorXd::Ones(x.size());
  p = coeffs(n)*ones;
  for (int j = n - 1; j >= 0; --j) {
    p = (x - t(j)*ones).cwiseProduct(p) + coeffs(j)*ones;
  }
}
