/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: December 2022
 */

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace VanDerPol {
/* SAM_LISTING_BEGIN_1 */
Eigen::Vector2d vanDerPolJacobianEVs(Eigen::Vector2d z, double mu) {
  Eigen::Matrix2d J(2, 2);
  J << 0, 1, -2 * mu * z[0] * z[1] - 1.0, mu * (1 - z[0] * z[0]);
  Eigen::EigenSolver<Eigen::Matrix2d> eig(J);
  auto ev = eig.eigenvalues().real();
  return ev.head(2);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
struct ROW {
  explicit ROW(unsigned int s)
      : s_(s),
        alpha_(Eigen::MatrixXd::Zero(s, s)),
        gamma_(Eigen::MatrixXd::Zero(s, s)),
        b_(Eigen::VectorXd::Zero(s)) {}
  unsigned int s_;         // Number of stages
  Eigen::MatrixXd alpha_;  // Coefficients $\cob{\alpha_{i,j}}$
  double d_gamma_;         // Diagonal of $\cob{\gamma}$-matrix
  Eigen::MatrixXd gamma_;  // Coefficients $\cob{\gamma_{i,j}}$
  Eigen::VectorXd b_;      // Weights $\cob{b_i}$
};
/* SAM_LISTING_END_2 */

// ROW coefficients for ROSxPR, J.Rang, J. Comp. Appl. Math. 286 (2015), 128-144
/* SAM_LISTING_BEGIN_3 */
ROW init_ROSxPR(void) {
  ROW row(3);                               // $\cob{s=3}$
  row.alpha_(1, 0) = 2.3660254037844388;    // $\cob{\alpha_{2,1}}$
  row.alpha_(2, 0) = 0.0;                   // $\cob{\alpha_{3,1}}$
  row.alpha_(2, 1) = 1.0;                   // $\cob{\alpha_{3,2}}$
  row.gamma_(1, 0) = -2.3660254037844388;   // $\cob{\gamma_{2,1}}$
  row.gamma_(2, 0) = -0.28468642516567449;  // $\cob{\gamma_{3,1}}$
  row.gamma_(2, 1) = -1.0813389786187642;   // $\cob{\gamma_{3,2}}$
  row.d_gamma_ = 0.78867513459481287;  // $\cob{\gamma_{1,1}=\gamma_{2,2}=\gamma_{3,3}}$
  row.b_[0] = 0.29266384402395124;    // $\cob{b_{1}}$
  row.b_[1] = -0.081338978618764143;  // $\cob{b_{2}}$
  row.b_[2] = 0.78867513459481287;    // $\cob{b_{3}}$
  return row;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
template <typename STATETYPE, typename RHSTYPE, typename JACTYPE>
std::vector<STATETYPE> odeROW(RHSTYPE &&rhs, JACTYPE &&Jacobian,
                              const STATETYPE &y0, double T, unsigned int M,
                              const ROW &row) {
  using MATTYPE = decltype(Jacobian(y0));  // Type for Jacobi matrix
  unsigned int N = y0.size();              // Dimension of state space
  std::vector<STATETYPE> y(M + 1);         // Vector of computed states
  y[0] = y0;
  std::vector<STATETYPE> k(row.s_);  // Increments $\cob{\Vk_i}$
  std::vector<STATETYPE> z(row.s_);  // Auxiliary states $\cob{\Vz_i}$
  const double h = T / M;            // timestep size
  // Main timestepping loop, y[m] contains the current iterate
  for (unsigned int m = 0; m < M; ++m) {
    const MATTYPE J = Jacobian(y[m]);
    const MATTYPE IgJ = MATTYPE::Identity(N, N) - h * row.d_gamma_ * J;
    auto IgJ_lu = IgJ.lu();  // Precompte LU decomposition
    z[0] = y[m];
    k[0] = IgJ_lu.solve(rhs(z[0]));          // Just fwd./bkw. substitution
    y[m + 1] = y[m] + h * row.b_[0] * k[0];  // Next state
    for (unsigned int i = 1; i < row.s_; ++i) {
      // Compute auxiliary states
      z[i] = y[m];
      for (unsigned int j = 0; j < i; ++j) {
        z[i] += h * row.alpha_(i, j) * k[j];
      }
      // Compute increments
      STATETYPE tmp = row.gamma_(i, 0) * k[0];
      for (unsigned int j = 1; j < i; ++j) {
        tmp += row.gamma_(i, j) * k[j];
      }
      k[i] =
          IgJ_lu.solve(rhs(z[i]) + h * J * tmp);  // Just fwd./bkw. substitution
      // Update next state
      y[m + 1] += h * row.b_[i] * k[i];
    }
  }
  return y;
}
/* SAM_LISTING_END_4 */

// Overload for scalar problems
/* SAM_LISTING_BEGIN_5 */
template <typename RHSTYPE, typename JACTYPE>
std::vector<double> odeROW(RHSTYPE &&rhs, JACTYPE &&Jacobian, double y0,
                           double T, unsigned int M, const ROW &row) {
  // Vector with computed states
  std::vector<double> y(M + 1);
  y[0] = y0;
  // Increments $\cob{\Vk_i}$
  std::vector<double> k(row.s_);
  // Auxiliary states $\cob{\Vz_i}$
  std::vector<double> z(row.s_);
  const double h = T / M;  // timestep size
  // Main timestepping loop, y[m] contains the current iterate
  for (unsigned int m = 0; m < M; ++m) {
    const double J = Jacobian(y[m]);
    const double IgJ = 1.0 - h * row.d_gamma_ * J;
    z[0] = y[m];
    k[0] = rhs(z[0]) / IgJ;
    y[m + 1] = y[m] + h * row.b_[0] * k[0];
    for (unsigned int i = 1; i < row.s_; ++i) {
      // Compute auxiliary states
      z[i] = y[m];
      for (unsigned int j = 0; j < i; ++j) {
        z[i] += h * row.alpha_(i, j) * k[j];
      }
      // Compute increments
      double tmp = row.gamma_(i, 0) * k[0];
      for (unsigned int j = 1; j < i; ++j) {
        tmp += row.gamma_(i, j) * k[j];
      }
      k[i] = (rhs(z[i]) + h * J * tmp) / IgJ;
      // Update next state
      y[m + 1] += h * row.b_[i] * k[i];
    }
  }
  return y;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
std::complex<double> stabFnROW(std::complex<double> z, const ROW &row) {
  std::vector<std::complex<double>> kv(row.s_);  // Increments $\cob{\Vk_i}$
  std::vector<std::complex<double>> zv(row.s_);  // States $\cob{\Vz_i}$
  const std::complex<double> IgJ = 1.0 - row.d_gamma_ * z;
  zv[0] = 1.0;
  kv[0] = z / IgJ;
  std::complex<double> S = 1.0 + row.b_[0] * kv[0];
  for (unsigned int i = 1; i < row.s_; ++i) {
    // Compute auxiliary states
    zv[i] = 1.0;
    for (unsigned int j = 0; j < i; ++j) {
      zv[i] += row.alpha_(i, j) * kv[j];
    }
    // Compute increments
    std::complex<double> tmp = row.gamma_(i, 0) * kv[0];
    for (unsigned int j = 1; j < i; ++j) {
      tmp += row.gamma_(i, j) * kv[j];
    }
    kv[i] = z * (zv[i] + tmp) / IgJ;
    S += row.b_[i] * kv[i];
  }
  return S;
}
/* SAM_LISTING_END_6 */

std::vector<double> rowDecay(double lambda, double T, unsigned int M) {
  // Initialize RWO method
  const ROW row = init_ROSxPR();
  // Define autonomous initial-value problem
  auto f = [lambda](double z) -> double { return lambda * z; };
  auto Df = [lambda](double /*z*/) -> double { return lambda; };
  double y0 = 1.0;
  // Launch timestepping
  return odeROW(f, Df, y0, T, M, row);
}

/* SAM_LISTING_BEGIN_7 */
std::vector<double> solveROWVanDerPol(double mu, double T, unsigned int M) {
  // Initialize ROW method
  const ROW row = init_ROSxPR();
  // Define autonomous initial-value problem
  auto f = [mu](Eigen::Vector2d z) -> Eigen::Vector2d {
    return Eigen::Vector2d(z[1], mu * (1.0 - z[0] * z[0]) * z[1] - z[0]);
  };
  auto Df = [mu](Eigen::Vector2d z) -> Eigen::Matrix2d {
    Eigen::Matrix2d J;
    J << 0, 1, -2 * mu * z[0] * z[1] - 1.0, mu * (1 - z[0] * z[0]);
    return J;
  };
  Eigen::Vector2d y0(2.0, 0.0);  // Initial state
  // Launch timestepping
  auto z_vec = odeROW(f, Df, y0, T, M, row);
  // The first component of the state vectors provides the y state
  std::vector<double> y_vec(z_vec.size());
  for (unsigned int k = 0; k < y_vec.size(); ++k) {
    y_vec[k] = z_vec[k][0];
  }
  return y_vec;
}
/* SAM_LISTING_END_7 */

std::vector<Eigen::Vector2d> rowVanDerPol(double mu, double T, unsigned int M) {
  // Initialize RWO method
  const ROW row = init_ROSxPR();
  // Define autonomous initial-value problem
  auto f = [mu](Eigen::Vector2d z) -> Eigen::Vector2d {
    return Eigen::Vector2d(z[1], mu * (1.0 - z[0] * z[0]) * z[1] - z[0]);
  };
  auto Df = [mu](Eigen::Vector2d z) -> Eigen::Matrix2d {
    Eigen::Matrix2d J;
    J << 0, 1, -2 * mu * z[0] * z[1] - 1.0, mu * (1 - z[0] * z[0]);
    return J;
  };
  Eigen::Vector2d y0(2.0, 0.0);
  // Launch timestepping
  return odeROW(f, Df, y0, T, M, row);
}

// Test convergence of single-step method based on error at final time
template <typename SSMEXP>
void testCVGSSM(SSMEXP &&solver, unsigned int M_min = 16,
                unsigned int no_runs = 8, unsigned int fac_ref = 10) {
  using STATETYPE = decltype(solver(1));
  // Vectors for storing final approximate states
  std::vector<STATETYPE> yT;
  unsigned int M = M_min;
  for (unsigned int l = 0; l < no_runs; ++l, M *= 2) {
    std::cout << "Solving with " << M << " timesteps" << std::endl;
    yT.push_back(solver(M));
  }
  // Finally run solver with very many timesteps to obtain reference solution.
  const STATETYPE y_ref = solver(fac_ref * M);
  // Tabulate the Euclidean vector norms of the final errors
  M = M_min;
  double y_errp;
  for (unsigned int l = 0; l < yT.size(); ++l, M *= 2) {
    const double y_errn = (y_ref - yT[l]).norm();
    std::cout << "M= " << M << ": |yM-y(T)| = " << y_errn;
    if (l > 0) {
      std::cout << ", ratio = " << y_errp / y_errn;
    }
    std::cout << std::endl;
    y_errp = y_errn;
  }
}

/* SAM_LISTING_BEGIN_9 */
void testCvgROWVanDerPol(void) {
  unsigned int M_min = 128;   // Minimal number of timesteps
  unsigned int no_runs = 8;   // Number of times the timestep is doubled 
  unsigned int fac_ref = 10;  // That many more timesteps for reference 
  // Vectors for storing final approximate states
  std::vector<double> yT;
  unsigned int M = M_min;
  for (unsigned int l = 0; l < no_runs; ++l, M *= 2) {
    yT.push_back(solveROWVanDerPol(1.0, 6.0, M).back());
  }
  // Finally run solver with very many timesteps to obtain reference solution.
  const double y_ref = solveROWVanDerPol(1.0, 6.0, fac_ref * M).back();
  // Tabulate the absolute values of the errors at final time
  M = M_min;
  double y_errp;
  for (unsigned int l = 0; l < yT.size(); ++l, M *= 2) {
    const double y_errn = std::abs(y_ref - yT[l]);
    std::cout << "M= " << M << ": |yM-y(T)| = " << y_errn;
    if (l > 0) {
      std::cout << ", ratio = " << y_errp / y_errn;
    }
    std::cout << std::endl;
    y_errp = y_errn;
  }
}
/* SAM_LISTING_END_9 */

}  // namespace VanDerPol

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "NumCSE code: Simulation of van der Pol equation" << std::endl;
  /* SAM_LISTING_BEGIN_X */
  std::cout << "Spectrum of Jacobian at [2;0] = "
            << VanDerPol::vanDerPolJacobianEVs(Eigen::Vector2d(2, 0), 1000)
                   .transpose()
            << std::endl;
  /* SAM_LISTING_END_X */

  std::cout << "Solution of Decay equation" << std::endl;
  auto states = VanDerPol::rowDecay(-1.0, 1.0, 20);
  for (auto y : states) {
    std::cout << y << ", ";
  }
  std::cout << std::endl;

  // Test stability function
  std::cout << "Test of stability function" << std::endl;
  std::vector<double> zv{1.0, -1.0, 2.0, -2.0, -3.0};
  for (double z : zv) {
    std::cout << "z= " << z
              << ": y_1 = " << VanDerPol::rowDecay(z, 1.0, 1).back()
              << " <-> S(z) = "
              << VanDerPol::stabFnROW(z, VanDerPol::init_ROSxPR()) << std::endl;
  }

  // Test solvers
  std::cout << "\n Simple decay equation" << std::endl;
  VanDerPol::testCVGSSM([](unsigned int M) -> Eigen::Matrix<double, 1, 1> {
    const double y_final = VanDerPol::rowDecay(-1.0, 1.0, M).back();
    return (Eigen::Matrix<double, 1, 1>() << y_final).finished();
  });

  // Test solvers
  std::cout << "\n Van der Pol equation" << std::endl;
  VanDerPol::testCVGSSM(
      [](unsigned int M) -> Eigen::Vector2d {
        return VanDerPol::rowVanDerPol(1.0, 10, M).back();
      },
      128);

  // Final test
  std::cout << "Empiric study of convergence" << std::endl;
  VanDerPol::testCvgROWVanDerPol();

  return 0;
}
