///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Till Ehrengruber <tille@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include "broyd.hpp"
#include "upbroyd.hpp"

int main() {
  // The function whose zero we want to find 
  auto F = [](const Eigen::Vector2d& x) -> Eigen::Vector2d {
    return Eigen::Vector2d(pow(x[0], 2) - pow(x[1], 4), x[0] - pow(x[1], 3));
  };
  // Jacobian of F 
  auto DF = [](const Eigen::Vector2d& x) -> Eigen::Matrix2d {
    Eigen::Matrix2d res;
    res << 2. * x[0], -4. * pow(x[1], 3), 1, -3. * pow(x[1], 2);
    return res;
  };
  // this function is called in every iteration and prints the progress
  auto print_progress_cb = [](unsigned k, Eigen::Vector2d x, Eigen::Vector2d f,
                              Eigen::Vector2d s) {
    std::cout << "Iteration " << k << ": |s| = " << s.norm()
              << ", |f(x)| =  " << f.norm() << std::endl;
  };
  // Initial guess
  Eigen::Vector2d x0(0.7, 0.7);
  {
    std::cout << "Matrix-based implementation of quasi-Newton method" << std::endl;
    const Eigen::Vector2d x = broyd(F, x0, DF(x0), 0.000001, 20, print_progress_cb);
    std::cout << "solution = " << x.transpose() << std::endl;
  }


  {
    // run the algorithm (we ignore the result)
    std::cout << "Recursion-based implementation of quasi-Newton method" << std::endl;
    const Eigen::Vector2d x = upbroyd(F, x0, DF(x0), 0.000001, 20, print_progress_cb);
    std::cout << "solution = " << x.transpose() << std::endl;
  }

  return 0;
}
