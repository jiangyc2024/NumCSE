/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

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

#include "helper_functions.hpp"

int main() {
  std::cout << "NumCSE Problem: periodic collocation equation \n" << std::endl;
  {
    constexpr unsigned int N = 10;
    constexpr unsigned int M = 16;
    const Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 1, 0, N);
    const Eigen::VectorXd y1 = eval_uN(x, M);
    const Eigen::VectorXd y2 = eval_uN_loop(x, M);
    std::cout << "eval_uN([" << x.transpose() << "], " << M
              << ") = " << y1.transpose() << std::endl
              << std::endl;
    std::cout << "Loop test: diff norm = " << (y1 - y2).norm() << std::endl
              << std::endl;
  }
  {
    constexpr unsigned int N = 10;
    const Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 1, 0, N);
    const Eigen::VectorXd f = eval_F(x);
    std::cout << "eval_F([" << x.transpose() << "]) = " << f.transpose()
              << std::endl
              << std::endl;
    std::cout << "eval_DF([" << x.transpose() << "]) = \n"
              << eval_DF(x) << std::endl;
  }

  return 0;
}
