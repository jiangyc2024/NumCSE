///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Till Ehrengruber <tille@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include "broyd.hpp"

int main() {
  {
    // The function whose zero we want to find
    auto F = [](const Eigen::Vector2d &x) -> Eigen::Vector2d {
      return Eigen::Vector2d(pow(x[0], 2) - pow(x[1], 4), x[0] - pow(x[1], 3));
    };
    // Jacobian of F
    auto DF = [](const Eigen::Vector2d &x) -> Eigen::Matrix2d {
      Eigen::Matrix2d res;
      res << 2. * x[0], -4. * pow(x[1], 3), 1, -3. * pow(x[1], 2);
      return res;
    };
    // this function is called in every iteration and prints the progress
    auto print_progress_cb = [](unsigned k, const Eigen::Vector2d & /*x*/,
                                const Eigen::Vector2d &f,
                                const Eigen::Vector2d &s) {
      std::cout << "Iteration " << k << ": |s| = " << s.norm()
                << ", |f(x)| =  " << f.norm() << std::endl;
    };
    // Initial guess
    Eigen::Vector2d x0(0.7, 0.7);
    {
      std::cout << "Matrix-based implementation of quasi-Newton method"
                << std::endl;
      const Eigen::Vector2d x =
          broyd(F, x0, DF(x0), 0.000001, 1.0E-9, 20, print_progress_cb);
      std::cout << "solution = " << x.transpose() << std::endl;
    }

    {
      std::cout << "Recursion-based implementation of quasi-Newton method"
                << std::endl;
      const Eigen::Vector2d x =
          upbroyd(F, x0, DF(x0), 0.000001, 1.0E-9, 20, print_progress_cb);
      std::cout << "solution = " << x.transpose() << std::endl;
    }
  }

  {
    std::cout
        << "TEST: Quasi-Newton method for large non-linear system of equations"
        << std::endl;
    constexpr unsigned int n = 100;
    const Eigen::VectorXd a{Eigen::VectorXd::LinSpaced(n, 0.0, n - 1) /
                            std::sqrt(0.5 * n * (n - 1) - 1.0)};
    const Eigen::MatrixXd A =
        Eigen::MatrixXd::Identity(n, n) + a * a.transpose();
    auto F_sub = [&A](const Eigen::VectorXd &x) -> Eigen::VectorXd {
      return (x.array() * ((A * x).array())).matrix();
    };
    // Manufactured solution 
    const Eigen::VectorXd sol{Eigen::VectorXd::Constant(n, 1.0)};
    const Eigen::VectorXd b = F_sub(sol);
    // Definition of F 
    auto F = [&F_sub, &b](const Eigen::VectorXd &x) -> Eigen::VectorXd {
      return (F_sub(x) - b);
    };
    // Definition of Jacobian of F
    auto DF = [&A](const Eigen::VectorXd &x) -> Eigen::MatrixXd {
      const auto xdA = x.asDiagonal() * A;
      const Eigen::MatrixXd Axd = (A * x).asDiagonal();
      return (xdA + Axd);
    };
    // Monitor for progress of iteration 
    auto monitor = [&sol](unsigned k, const Eigen::VectorXd &x,
                          const Eigen::VectorXd &f, const Eigen::VectorXd &s) {
      std::cout << "Quasi-Newton iteration " << k << ": |s| = " << s.norm()
                << ", |f(x)| =  " << f.norm()
                << ", |err| = " << (x - sol).norm() << std::endl;
    };
    // Initial guess
    const Eigen::VectorXd x0{Eigen::VectorXd::LinSpaced(n, 0.5, 1.0)};

    constexpr double reltol = 1.0E-6;
    constexpr double abstol = 1.0E-9;
    std::cout << "Quasi-Newton for " << x0.size() << " x " << x0.size()
              << " non-linear system of equations" << std::endl;
    std::cout << "Non-recursive implementation" << std::endl;
    broyd(F, x0, DF(x0), reltol, abstol, 20, monitor);
    std::cout << "Recursive implementation" << std::endl;
    upbroyd(F, x0, DF(x0), reltol, abstol, 20, monitor);
  }

  return 0;
}
