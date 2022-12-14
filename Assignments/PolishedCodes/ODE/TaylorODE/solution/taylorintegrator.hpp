#ifndef TAYLORINTEGRATORHPP
#define TAYLORINTEGRATORHPP

#include <Eigen/Dense>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include "ode45.hpp"

/**
 * \brief Solves the Predator Prey model using the Taylor expansion method.
 *
 * \param alpha1 parameter of ODE
 * \param beta1 parameter of ODE
 * \param alpha2 parameter of ODE
 * \param beta2 parameter of ODE
 * \param T final time
 * \param y0 initial state
 * \param M number of steps to take
 * \return std::vector<Eigen::Vector2d> of states
 */
/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Vector2d> SolvePredPreyTaylor(double alpha1, double beta1,
                                                 double alpha2, double beta2,
                                                 double T,
                                                 const Eigen::Vector2d &y0,
                                                 unsigned int M) {
  // Vector storing the states
  std::vector<Eigen::Vector2d> res;
  // TODO: (11-8.c) Solve the predator prey model using the Taylor expansion
  // method.
  // START
  res.reserve(M + 1);

  // Step size
  const double h = T / M;

  // Storing the initial state
  Eigen::Vector2d y = y0;
  res.push_back(y);

  // Preparing functions for using the Taylor expansion method
  // Lambda function for the RHS f(y)
  auto f = [&](const Eigen::Vector2d &y) {
    Eigen::Vector2d f;
    f << (alpha1 - beta1 * y(1)) * y(0), (beta2 * y(0) - alpha2) * y(1);
    return f;
  };

  // Lambda function for computing Df(y)*z
  auto Df = [&](const Eigen::Vector2d &y, const Eigen::Vector2d &z) {
    Eigen::Matrix2d Df;
    Df << alpha1 - beta1 * y(1), -beta1 * y(0), beta2 * y(1),
        -alpha2 + beta2 * y(0);
    Eigen::Vector2d out = Df * z;
    return out;
  };

  // Lambda function for D^2f(y) (z,z)
  auto D2f = [&](const Eigen::Vector2d &y, const Eigen::Vector2d &z) {
    Eigen::Matrix2d H1;
    Eigen::Matrix2d H2;
    H1 << 0, -beta1, -beta1, 0;
    H2 << 0, beta2, beta2, 0;
    Eigen::Vector2d out;
    out << z.transpose() * H1 * z, z.transpose() * H2 * z;
    return out;
  };

  for (unsigned int k = 0; k < M; ++k) {
    // evaluate terms for taylor step.
    auto fy = f(y);
    // Df(y) * f(y)
    auto dfyfy = Df(y, fy);
    // Df(y) * Df(y) * f(y)
    auto df2yfy = Df(y, dfyfy);
    // D^2f(y) (f(y),f(y))
    auto d2fyfy = D2f(y, fy);

    // evaluate taylor expansion to compute update
    y = y + h * fy + 0.5 * h * h * dfyfy +
        1.0 / 6.0 * h * h * h * (df2yfy + d2fyfy);

    // save new state:
    res.push_back(y);
  }
  // END
  return res;
}
/* SAM_LISTING_END_1 */

void PrintErrorTable(const Eigen::ArrayXd &M, const Eigen::ArrayXd &error) {
  std::cout << std::setw(15) << "M" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;

  for (unsigned int i = 0; i < M.size(); ++i) {
    std::cout << std::setw(15) << M(i) << std::setw(15) << error(i);
    if (i > 0) {
      std::cout << std::setw(15) << std::log2(error(i - 1) / error(i));
    }
    std::cout << std::endl;
  }
}

/* SAM_LISTING_BEGIN_2 */
double testCvgTaylorMethod() {
  double cvgrate = 0;
  // TODO: (11-8.d) Generate an error table and return an estimate of the
  // convergence rate of the Taylor method.
  // START
  // initialize parameters for the model:
  constexpr double T = 10;     // final time
  Eigen::Vector2d y0(100, 5);  // initial condition
  Eigen::Vector2d yex(0.319465882659820,
                      9.730809352326228);  // reference solution
  constexpr double alpha1 = 3.0;
  constexpr double alpha2 = 2.0;
  constexpr double beta1 = 0.1;
  constexpr double beta2 = 0.1;

  // Initialize parameters for the convergence study
  constexpr unsigned int M0 = 128;    // Minimum number of timesteps
  constexpr unsigned int numRef = 8;  // Number of refinements

  // Convergence study
  Eigen::ArrayXd error(numRef);
  Eigen::ArrayXd M(numRef);
  for (unsigned int i = 0; i < numRef; ++i) {
    M(i) = std::pow(2, i) * M0;
    auto res = SolvePredPreyTaylor(alpha1, beta1, alpha2, beta2, T, y0, M(i));
    error(i) = (res.back() - yex).norm();
  }

  PrintErrorTable(M, error);

  // calculate linear regression line: log(error) ~ c0 + c1*log(M)
  Eigen::MatrixXd A(numRef, 2);
  A.col(0) = Eigen::VectorXd::Ones(numRef);
  A.col(1) = M.log();
  Eigen::VectorXd logError = error.log();
  Eigen::Vector2d coeffs = A.householderQr().solve(logError);

  // estimated convergence rate: -c1
  cvgrate = -coeffs(1);
  // END
  return cvgrate;
}
/* SAM_LISTING_END_2 */

#endif
