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

template <typename Functor, typename Jacobian>
std::array<Eigen::VectorXd, 5>
computeStages(Functor &&f, Jacobian &&J, const Eigen::VectorXd &y, double h,
              double rtol = 1E-6, double atol = 1E-8) {
  std::array<Eigen::VectorXd, 5> G; // array of stages 
  int d = y.size();
  MatrixXd Coeffs = ButcherMatrix();
  // TO DO (0-2.d)
  // START
  
  // END
  return G;
}

template <typename Functor, typename Jacobian>
Eigen::VectorXd discEvolSDIRK(Functor &&f, Jacobian &&J,
                              const Eigen::VectorXd &y, double h,
                              double rtol = 1E-6, double atol = 1E-8) {
  // The b weights are in the last row of Coeffs.
  MatrixXd Coeffs = ButcherMatrix();
  VectorXd Psi;
  // TO DO (0-2.e)
  // START
  
  // END
  return Psi;
}

std::vector<Eigen::VectorXd> solveGradientFlow(const Eigen::VectorXd &d,
                                               double lambda,
                                               const Eigen::VectorXd &y,
                                               double T, unsigned int N) {
  std::vector<Eigen::VectorXd> Y(N+1);
  // TO DO (0-2.h)
  // START
  
  // END
  return Y;
}

