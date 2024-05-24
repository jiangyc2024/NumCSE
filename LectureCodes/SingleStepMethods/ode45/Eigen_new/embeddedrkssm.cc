/**
 * @file
 * @brief function realizing an embedded Runge-Kutta-Fehlberg explicit
 * single-step method
 * @author Ralf Hiptmair
 * @date   April 2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

namespace EmbeddedRKSSM {

struct Options {
  double T;       // final time
  double h0;      // initial stepsize
  double reltol;  // reltol relative tolerance for timestep control
  double abstol;  // abstol absolute tolerance for timestep control
  double hmin;    // minimal admissible stepsize
} __attribute__((aligned(64)));

/**
 * @brief Adaptive embedded Runge-Kutta-Fehlberg timestepping
 *
 * @tparam RHSFunction functor type for right-hand-side vector field
 * @param f_rhs functor providing the right-hand-side function of the autonomous
 * ODE through its evaluation operator
 * @param A Butcher matrix of size s x s (only strictly lower triangular part
 * used)
 * @param b weight vector for RKSSM of order p+1
 * @param bh weight vector for RKSSM of order p
 * @param p order of lower-order method
 * @param y0 intial value
 * @param opt user options
 */
template <typename RHSFunction>
std::vector<std::pair<double, Eigen::VectorXd>> embeddedRKSSM(
    RHSFunction &&f_rhs, const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
    const Eigen::VectorXd &bh, unsigned int p, const Eigen::VectorXd &y0,
    const Options &opt) {
  // Check parameters defining embedded RK-SSM
  const unsigned int s = A.cols();  // Number of stages
  assert((s == A.rows()) && "Butcher matrix must be square");
  assert((s == b.size()) && "Length of weight vector b != no of stages");
  assert((s == bh.size()) && "Length of weight vector bh != no of stages");
  const unsigned int N = y0.size();  // Dimension of state space
  Eigen::MatrixXd K(N, s);     // Columns hold increment vectors

  double t = 0.0;  // Initial time zero for autonomous initial-value problem
  double h = opt.h0;   // Current timestep size
  std::vector<std::pair<double, Eigen::VectorXd>> states{
      {t, y0}};            // State sequence
  Eigen::VectorXd y = y0;  // Current state
  states.emplace_back(t, y);
  // Main timestepping loop
  while ((states.back().first < opt.T) && (h >= opt.hmin)) {
    // Compute increments
    K.col(0) = f_rhs(y);
    for (unsigned int l = 1; l < s; ++l) {
      Eigen::VectorXd v{Eigen::VectorXd::Zero(N)};
      for (unsigned int i = 0; i < l; ++i) {
        v += A(l, i) * K.col(i);
      }
      K.col(l) = f_rhs(y + h * v);
    }
    // Compute next two approximate states
    auto yh = y + h * K * b;   // high-order method
    auto yH = y + h * K * bh;  // low-order method
    const double est = (yh - yH).norm();
    const double tol = std::max(opt.reltol * y.norm(), opt.abstol);
    if (est <= tol) {
      y = yh;  // Advance to next approximate state
      t += h;  // Next time
      states.emplace_back(t, y);
      // std::cout << "t = " << t << ", y = " << y.transpose() << std::endl;
    }
    // else {
    //   std::cout << "tol/est = " << tol / est << ", h = " << h << std::endl;
    // }
    h *= std::max(0.5, std::min(2., 0.9 * std::pow(tol / est, 1. / (p + 1))));
    if (h < opt.hmin) {
      std::cerr
          << "Warning: Failure at t=" << states.back().first
          << ". Unable to meet integration tolerances without reducing the step"
          << " size below the smallest value allowed (" << opt.hmin
          << ") at time t." << std::endl;
    } else {
      h = std::min(opt.T - t + opt.hmin, h);
    }
  }
  return states;
}

// Helper function for testing
void testrun(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
             const Eigen::VectorXd &bh, unsigned int p) {
  std::cout << "Test run of embedded RK-SSM" << std::endl;
  std::cout << "Butcher matrix = \n " << A << std::endl;
  std::cout << "Weight vector (order p+1) = " << b.transpose() << std::endl;
  std::cout << "Weight vector (order p) = " << bh.transpose() << std::endl;
  // Simple linear test case: rotation ODE
  auto f = [](Eigen::Vector2d y) -> Eigen::Vector2d {
    return {-y[1], y[0]};
  };
  Eigen::VectorXd y0(2);
  y0 << 1.0, 0.0;
  // Final time
  const double T = 2.0 * 3.14159265358979323846;
  // Test different tolerances
  const int m = 6;
  std::vector<double> rtol{0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};
  std::vector<double> atol{0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001};
  for (int j = 0; j < m; ++j) {
    auto res = embeddedRKSSM(f, A, b, bh, p, y0, {.T=T, .h0= T / 100, .reltol=rtol[j], .abstol=atol[j], .hmin=T / 1E6});
    double err = 0.0;
    for (const auto &i : res) {
      const double t = i.first;
      const Eigen::Vector2d exact(std::cos(t), std::sin(t));
      const Eigen::Vector2d approx = i.second;
      err = std::max(err, (exact - approx).norm());
    }
    std::cout << "rtol = " << rtol[j] << ", atol = " << atol[j] << " : "
              << (res.size() - 1) << " steps, err = " << err << std::endl;
  }
}

}  // namespace EmbeddedRKSSM

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "Adaptive embedded Runge-Kutta-Fehlberg single-step method"
            << std::endl;
  {
    // Simplest embedded method: Euler - Heun
    Eigen::MatrixXd A(2, 2);
    A << 0, 0, 1, 0;
    Eigen::VectorXd b(2);
    b << 0.5, 0.5;
    Eigen::VectorXd bh(2);
    bh << 1.0, 0.0;
    EmbeddedRKSSM::testrun(A, b, bh, 1);
  }

  {
    // More complicated: Bogacki-Shampine
    // https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Embedded_methods
    Eigen::MatrixXd A(4, 4);
    A << 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.75, 0, 0, 2. / 9., 1. / 3., 4. / 9., 0;
    Eigen::VectorXd b(4);
    b << 2.0 / 9., 1. / 3., 4. / 9., 0;
    Eigen::VectorXd bh(4);
    bh << 7. / 24., 1. / 4., 1. / 3., 1. / 8.;
    EmbeddedRKSSM::testrun(A, b, bh, 2);
  }
  return 0;
}
