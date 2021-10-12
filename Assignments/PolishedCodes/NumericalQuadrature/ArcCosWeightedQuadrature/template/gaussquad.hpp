/* **********************************************************************
 * Course "Numerical Methods for CSE", R. Hiptmair, SAM, ETH Zurich
 * Author: R. Hiptmair
 * Date: January 2020
 */

#ifndef GAUSSQUAD_HPP
#define GAUSSQUAD_HPP
#include <Eigen/Dense>

// A simple data structure for storing a quadrature rule
struct QuadRule {
  QuadRule() = default;
  explicit QuadRule(unsigned int n) : nodes_(n), weights_(n) {}
  Eigen::VectorXd nodes_;
  Eigen::VectorXd weights_;
};

QuadRule gaussquad(unsigned int n);
#endif
