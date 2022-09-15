/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

#ifndef HELPERFUNCTIONSHPP
#define HELPERFUNCTIONSHPP

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <unsupported/Eigen/FFT>
#include <vector>
#include "periodiccollocation.hpp"

// Generic Newton's method with correction-based termination based on simplified
// Newton correction
template <typename FuncType, typename JacType, typename VecType,
          typename RECORDER = std::function<void(const VecType &, double)>>
void newton_stc(
    const FuncType &F, const JacType &DF, VecType &x, double rtol, double atol,
    RECORDER &&rec = [](const VecType &, double) -> void {}) {
  std::cout << "NEWTON: Initial guess = " << x.transpose() << std::endl;
  using scalar_t = typename VecType::Scalar;
  scalar_t sn;
  do {
    auto jacfac = DF(x).lu();  // LU-factorize Jacobian \Label[line]{nstc:1}]
    x -= jacfac.solve(F(x));   // Compute next iterate \Label[line]{nstc:2}
    // Compute norm of simplified Newton correction
    sn = jacfac.solve(F(x)).norm();
    std::cout << "NEWTON: step, x = " << x.transpose() << std::endl;
    rec(x, sn);  // record iteration history
  }
  // Termination based on simplified Newton correction
  while ((sn > rtol * x.norm()) && (sn > atol));
}

Eigen::VectorXd eval_F_zero_rhs(const Eigen::VectorXd &x) {
  const unsigned int N = x.size() - 1;
  Eigen::VectorXd Fx;
  Eigen::ArrayXd s = Eigen::ArrayXd::LinSpaced(N + 1, 0, N);
  s = 4 * M_PI * M_PI * s * s;
  const Eigen::VectorXd d{(x.array() * s).matrix()};
  Fx = eval_uN(d, N + 1) + (eval_uN(x, N + 1).array().pow(3)).matrix();
  return Fx;
}

// Newton's method for non-linear system of equations arising from periodic
// collocation with N+1 cosine functions
Eigen::VectorXd solve_per_coll(unsigned int N) {
  // Lambda function realizing F
  auto F = [&N](Eigen::VectorXd &x) -> Eigen::VectorXd {
    // START student code
    return eval_F(x);
    // END student code
  };
  // Lambda function for the Jacobian
  auto DF = [&N](Eigen::VectorXd &x) -> Eigen::MatrixXd {
    return eval_DF(x);
  };
  // Initial guess (a constant function)
  Eigen::VectorXd x{Eigen::VectorXd::Constant(N + 1, 1.0)};
  // x[0] = 1.0;
  // Variables for recording progress of iteration
  std::vector<Eigen::VectorXd> rec_x;
  std::vector<double> rec_sn;
  auto rec = [&rec_x, &rec_sn](const Eigen::VectorXd &x, double sn) -> void {
    rec_x.push_back(x);
    rec_sn.push_back(sn);
  };
  newton_stc(F, DF, x, 1.0E-8 /*rtol*/, 1.0E-10 /*atol*/, rec);
  std::cout << std::setw(15) << " |x-x(k)| " << std::setw(15) << " sn "
            << std::endl;
  const int n_steps = rec_x.size();
  for (int l = 0; l < n_steps; ++l) {
    std::cout << std::setw(15) << std::scientific << (x - rec_x[l]).norm()
              << std::setw(15) << std::scientific << rec_sn[l] << std::endl;
  }
  return x;
}

// Loop based evaluation of trigonometric polynomial
Eigen::VectorXd eval_uN_loop(const Eigen::VectorXd &x, unsigned int M) {
  unsigned int N = x.size() - 1;
  Eigen::VectorXd y{Eigen::VectorXd::Zero(M)};
  for (unsigned int k = 0; k < M; ++k) {
    double tk = (double)k / (double)M;
    for (unsigned int j = 0; j <= N; ++j) {
      y[k] += x[j] * std::cos(2.0 * M_PI * j * tk);
    }
  }
  return y;
}

#endif
