#include <Eigen/Dense>
#include <array>
#include <iostream>
#include <vector>

using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd ButcherMatrix() {
  Eigen::MatrixXd A(6, 5);
  // clang-format off
  A <<       0.25,          0.,       0.,       0.,   0.,
    0.5,        0.25,       0.,       0.,   0.,
    17./50.,     -1./25.,     0.25,       0.,   0.,
    371./1360., -137./2720., 15./544.,     0.25,   0.,
    25./24.,    -49./48., 125./16., -85./12., 0.25,
    25./24.,    -49./48., 125./16., -85./12., 0.25;
  // clang-format on
  return A;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <typename Functor, typename Jacobian>
Eigen::VectorXd solveGenStageEquation(Functor &&f, Jacobian &&J,
                                      const Eigen::VectorXd &y,
                                      const Eigen::VectorXd &b, double h,
                                      double rtol = 1E-6, double atol = 1E-8) {
  // Need to solve the equation lhs(g) = g - h*f(y+g)/4 - b = 0.
  auto lhs = [f, y, b, h](const VectorXd &g) {
    // lhs(g) = g - h*f(y+g)/4 - b
    VectorXd val = g - 0.25 * h * f(y + g) - b;
    return val;
  };
  auto Jlhs = [J, y, h](const VectorXd &g) {
    // lhs(g) = g - h*f(y+g)/4 - b, so the Jacobian is
    // Jlhs(g) = Id - h*Jf(y+g)/4
    int dim = y.size();
    MatrixXd Jval = MatrixXd::Identity(dim, dim) - 0.25 * h * J(y + g);
    return Jval;
  };

  // Perform Newton iterations:
  VectorXd g = VectorXd::Zero(y.size());        // initial guess g=0.
  VectorXd delta = -Jlhs(g).lu().solve(lhs(g)); // Newton correction term.
  int iter = 0, maxiter = 100; // If correction based termination does not work.
  while (delta.norm() > atol && delta.norm() > rtol * g.norm() &&
         iter < maxiter) {
    g = g + delta;
    delta = -Jlhs(g).lu().solve(lhs(g));
    iter++;
  }
  return g + delta; // Perform the final step
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename Functor, typename Jacobian>
std::array<Eigen::VectorXd, 5> computeStages(Functor &&f, Jacobian &&J,
                                             const Eigen::VectorXd &y, double h,
                                             double rtol = 1E-6,
                                             double atol = 1E-8) {
  std::array<Eigen::VectorXd, 5> G; // array of stages
  int d = y.size();
  MatrixXd Coeffs = ButcherMatrix();
  // TO DO (12-8.e)
  // START

  // END
  return G;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_5 */
template <typename Functor, typename Jacobian>
Eigen::VectorXd discEvolSDIRK(Functor &&f, Jacobian &&J,
                              const Eigen::VectorXd &y, double h,
                              double rtol = 1E-6, double atol = 1E-8) {
  // The b weights are in the last row of Coeffs.
  MatrixXd Coeffs = ButcherMatrix();
  VectorXd Psi;
  // TO DO (12-8.f)
  // START

  // END
  return Psi;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
std::vector<Eigen::VectorXd> solveGradientFlow(const Eigen::VectorXd &d,
                                               double lambda,
                                               const Eigen::VectorXd &y,
                                               double T, unsigned int M) {
  std::vector<Eigen::VectorXd> Y(M+1);
  // TO DO (12-8.i)
  // START
  
  // END
  return Y;
}
/* SAM_LISTING_END_6 */
