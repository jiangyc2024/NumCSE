/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */
 
#ifdef NICEBACKTRACE
#include "backtrace.hpp"
#endif

#include "locationestimation.hpp"
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#ifdef NICEBACKTRACE
#include "backtrace.hpp"
#endif

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "NumCSE code for location estimation of an acoustic source\n"
            << std::endl;
  // Creation of synthetic test data
  // Initialize receiver positions at 5 corners of unit cube
  Eigen::Matrix<double, 3, Eigen::Dynamic> Q(3, 11);
  // clang-format off
  Q << 0, 0, 0, 1, 1, 1, 1, 0.5, 0.0, 1.0, 0.5,
       1, 0, 1, 1, 0, 0, 1, 0.0, 0.5, 0.5, 1.0,
       0, 1, 1, 0, 1, 0, 1, 0.5, 0.5, 0.5, 0.5;
  // clang-format on
  // Source point
  const Eigen::Vector3d ps{(Eigen::Vector3d() << 0.1, 0.1, 0.9).finished()};
  // Starting time
  double ts = 0.0;
  // Compute arrival times
  const int n = Q.cols();
  Eigen::VectorXd ta(n);
  for (int j = 0; j < n; ++j) {
    ta[j] = ts + (Q.col(j) - ps).norm();
  }
  // Estimate source location and starting time using Gauss-Newton Method
  std::vector<double> rec_err;
  std::vector<double> rec_sn;
  std::vector<double> rec_Fn;
  Eigen::Vector4d xs{(Eigen::Vector4d() << ps, ts).finished()};
  auto rec = [&rec_err, &rec_sn, &rec_Fn, &xs](const Eigen::Vector4d &x,
                                               const Eigen::Vector4d &s,
                                               double Fn) -> void {
    rec_err.push_back((x - xs).norm());
    rec_sn.push_back(s.norm());
    rec_Fn.push_back(Fn);
    std::cout << "GN step: |s| = " << s.norm() << ", |F(x)| = " << Fn
              << std::endl;
  };
  auto [p_est, t_est] = source_estimation(Q, ta, rec);
  std::cout << "\nEstimated location = " << p_est.transpose()
            << ", estimated time = " << t_est << std::endl;
  std::cout << "\nIteration history for Gauss-Newton method" << std::endl;
  std::cout << std::setw(15) << "k" << std::setw(15) << "error norm "
            << std::setw(15) << "|s|" << std::setw(15) << "|F(x)|" << std::endl;
  for (unsigned int l = 0; l < rec_err.size(); ++l) {
    std::cout << std::setw(15) << l + 1 << " & " << std::setw(15)
              << std::scientific << rec_err[l] << " & " << std::setw(15)
              << std::scientific << rec_sn[l] << " & " << std::setw(15)
              << std::scientific << rec_Fn[l] << "\\\\" << std::endl;
  }
  return 0;
}
