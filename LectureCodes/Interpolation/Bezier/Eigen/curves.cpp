///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2020 SAM, D-MATH
/// Author(s): R. Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <cassert>
#include <iostream>

/** @brief Evalutation of 2D Bezier polynomial in multiple points
    @param nodes 2xn matrix of control point coordinates
    @param t row vector of N parameter values for which to evaluate
    @return row vector of length N

     The function evaluates a Bezier polynomial curve defined through
     its control points for all parameters given in a vector.

     Algorithm: Modified Horner scheme, not De Casteljau

     @note Assumes the parameter interval [0,1]
*/
/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd evalBezier(const Eigen::MatrixXd &nodes,
                           const Eigen::RowVectorXd &t) {
  assert((nodes.rows() == 2) && "Nodes must have two coordinates");
  const Eigen::Index n = nodes.cols();  // Number of control points
  const Eigen::Index d = n - 1;         // Polynomial degree
  const Eigen::Index N = t.size(); // No. of evaluation points
  // Vector containing 1-t ("one minus t")
  const auto oml{Eigen::RowVectorXd::Constant(N, 1.0) - t};
  // Modified Horner scheme for polynomial in Bezier form
  // Vector for returning result, initialized with p[0]*(1-t)
  Eigen::MatrixXd res = nodes.col(0) * oml;
  double binom_val = 1.0; // Accumulate binomial coefficients
  // Powers of argument values
  Eigen::RowVectorXd t_pow{Eigen::RowVectorXd::Constant(N, 1.0)};
  for (int i = 1; i < d; ++i) {
    t_pow.array() *= t.array();
    binom_val *= (static_cast<double>(d - i) + 1.0) / i;
    res += binom_val * nodes.col(i) * t_pow;
    res.row(0).array() *= oml.array();
    res.row(1).array() *= oml.array();
  }
  res += nodes.col(d) * (t.array() * t_pow.array()).matrix();
  return res;
}
/* SAM_LISTING_END_1 */

int main(int /*argc*/, char ** /*argv*/) {
  // Initialize nodes
  // clang-format off
  const Eigen::MatrixXd nodes = ((Eigen::MatrixXd(2,9) <<
    -1. ,  0. ,  4. ,  5.5,  9.5, 11. , 15. , 15.5, 15.5,
     0. ,  2. ,  2. ,  4. ,  4. ,  2. ,  2. ,  3. ,  0.).finished());
  // clang-format on 
  const Eigen::RowVectorXd lambda = (Eigen::ArrayXd::LinSpaced(10,0.0,1.0)).matrix();
  const Eigen::MatrixXd curvevals = evalBezier(nodes,lambda);
  std::cout << "Points on curve" << std::endl << curvevals << std::endl;
  return 0;
}

// clang-format off
/* 
Expected result:

[-1.        ,  0.59801478,  2.83632398,  5.22473708,  7.63447843, 9.99453481, 12.19067215, 14.01608544, 15.17578928, 15.5       ],
[ 0.        ,  1.31909067,  2.20708583,  2.8117665 ,  3.0347603 ,  2.890281  ,  2.55570797,  2.19174808,  1.62490908,  0.        ]
*/
