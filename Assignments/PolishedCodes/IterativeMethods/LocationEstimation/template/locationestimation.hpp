/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2022
 */

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

#define _USE_MATH_DEFINES

// Function for solving the location estimation problem.
// The positions of the receivers are passed in the matrix P, the arrival times
// in the vector ta. Speed of sound is set to 1
template <typename RECORDER = std::function<
              void(const Eigen::Vector4d &, const Eigen::Vector4d &, double)>>
std::pair<Eigen::Vector3d, double> source_estimation(
    const Eigen::Matrix<double, 3, Eigen::Dynamic> &Q,
    const Eigen::VectorXd &ta,
    RECORDER &&rec = [](const Eigen::Vector4d &, const Eigen::Vector4d &,
                        double) -> void {}) {
  const int n = ta.size();
  assert(Q.cols() == n);
  // Parameters
  Eigen::Vector3d p;
  double t;
  // TODO : Write a function that solves the non linear system of equations 
  // (eq. 0.1.1 in exam PDF) in the least squares sense using Gauss Newton method
  // START
  
  // END
  return {p, t};
}
