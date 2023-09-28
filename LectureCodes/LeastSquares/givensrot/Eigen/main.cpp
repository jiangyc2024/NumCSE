///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>

#include <Eigen/Dense>

#include "givenscoltrf.hpp"
#include "planerot.hpp"
#include "qrgivens.hpp"

using Eigen::VectorXd;
using Eigen::MatrixXd;

int main() {
  std::cout << "#####################################" << std::endl;
  std::cout << "New test of planerot" << std::endl;
  std::cout << "#####################################" << std::endl;
  Eigen::Vector2d a;
  a << 2, 1;
  Eigen::Matrix2d G;
  Eigen::Vector2d x;
  planerot::planerot(a, G, x);
  std::cout << "a=\n"
            << a << std::endl
            << "G=\n"
            << G << std::endl
            << "G*G^T = \n"
            << G.transpose() * G << std::endl
            << "x=\n"
            << x << std::endl;
  std::cout << "Error norm (Gx-a) = " << (G.transpose() * a - x).norm()
            << std::endl;

  std::cout << "#####################################" << std::endl;
  std::cout << "Test of givenscoltrf" << std::endl;
  std::cout << "#####################################" << std::endl;
  const int n = 4;
  VectorXd r1(n);
  Eigen::VectorXd aIn(n);
  aIn << 1, 3, 4, 8; //, 9, 42, 2401, 343;
  MatrixXd Q1(n, n);
  givenscoltrf::givenscoltrf(aIn, Q1, r1);
  std::cout << "aIn=\n"
            << aIn << std::endl
            << "Q=\n"
            << Q1 << std::endl
            << "aOut=\n"
            << r1 << std::endl;
  std::cout << "Error norm (Q'Q-I) = "
            << (Q1.transpose() * Q1 - Eigen::MatrixXd::Identity(n, n)).norm()
            << std::endl;
  std::cout << "Error norm (Qr-a) = " << (Q1 * r1 - aIn).norm() << std::endl;

  std::cout << "#####################################" << std::endl;
  std::cout << "Test of qrgivens" << std::endl;
  std::cout << "#####################################" << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  MatrixXd AIn(n, n);
  AIn << 1, 2, 4, 7, 9, 4, 2, 9, 7, 5, 9, 2, 3, 4, 4, 8;
  MatrixXd R2(n, n);
  MatrixXd Q2(n, n);
  qrgivens::qrgivens(AIn, Q2, R2);
  std::cout << "AIn=\n"
            << AIn << std::endl
            << "Q=\n"
            << Q2 << std::endl
            << "R=\n"
            << R2 << std::endl;
  std::cout << "Error norm (Q'Q-I) = "
            << (Q2.transpose() * Q2 - Eigen::MatrixXd::Identity(n, n)).norm()
            << std::endl;
  std::cout << "Error norm (QR-A) = " << (Q2 * R2 - AIn).norm() << std::endl;
  return 0;
}
