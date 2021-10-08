////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "linfit.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main() {
  // Declare test data $(t_i, f_i)$
  Eigen::VectorXd t_vec = Eigen::VectorXd::LinSpaced(10, 0.1, 1.0);
  Eigen::VectorXd f(10);
  f << 100.0, 34.0, 17.0, 12.0, 9.0, 6.0, 5.0, 4.0, 4.0, 2.0;

  // Create matrix A
  Eigen::MatrixXd A = make_A(t_vec);

  std::cout << "\nThe matrix A is defined by \n\n" << A << std::endl << std::endl;

  // Compute coefficients using the normal equation
  Eigen::VectorXd gamma1 = data_fit_normal(t_vec, f);

  // Compute coefficients using QR decomposition
  Eigen::VectorXd gamma2 = data_fit_qr(t_vec, f);

  // Create time vector and function value vector used to plot the fitted
  // function.
  Eigen::VectorXd tl = Eigen::VectorXd::LinSpaced(91, 0.1, 1.0);
  Eigen::VectorXd yl1, yl2;

  Eigen::VectorXd err1 = fitting_error(gamma1);
  Eigen::VectorXd err2 = fitting_error(gamma2);

  // Approximate function coefficients using both methods
  std::cout << "Enter \"1\" to plot data_fit_normal().\n"
            << "Enter \"2\" to plot data_fit_qr().\n"
            << "Enter another key to exit.\n";
  int ans = 0;
  std::cin >> ans;

  switch (ans) {
    case 1:
      std::cout << "\nUsing data_fit_normal: gamma =\n" << gamma1 << std::endl;

      // Get function values. Needed to plot the fitted function.
      yl1 = fitted_function(gamma1, tl);

      // Plot data points and fitted function
      plt::figure();
      plt::semilogy(tl, yl1, {{"label", "normal equation"}});
      plt::semilogy(t_vec, f, "o", {{"label", "data"}});
      plt::title("Fitted function");
      plt::xlabel("t");
      plt::ylabel("y");
      plt::legend("best");
      plt::savefig("./cx_out/fit_normal.png");

      // Plot fitting errors
      plt::figure();
      plt::semilogy(t_vec, err1, "o", {{"label", "normal equation"}});
      plt::title("Fit error squared");
      plt::xlabel("t");
      plt::ylabel("|y-f(t)|^2");
      plt::legend("best");
      plt::savefig("./cx_out/error_normal.png");

      std::cout << "\nSee Files for plots\n";

      break;
    case 2:
      std::cout << "\nUsing data_fit_qr: gamma =\n" << gamma2 << std::endl;

      // Get function values. Needed to plot the fitted function.
      yl2 = fitted_function(gamma2, tl);

      // Plot data points and fitted function
      plt::figure();
      plt::semilogy(tl, yl2, {{"label", "qr fitting"}});
      plt::semilogy(t_vec, f, "o", {{"label", "data"}});
      plt::title("Fitted function");
      plt::xlabel("t");
      plt::ylabel("y");
      plt::legend("best");
      plt::savefig("./cx_out/fit_qr.png");

      // Plot fitting errors
      plt::figure();
      plt::semilogy(t_vec, err2, "o", {{"label", "qr fitting"}});
      plt::title("Fit error squared");
      plt::xlabel("t");
      plt::ylabel("|y-f(t)|^2");
      plt::legend("best");
      plt::savefig("./cx_out/error_qr.png");
      
      std::cout << "\nSee Files for plots\n";
      
      break;
    default:
      return 0;
  }

  // Show that the different methods don't give the same result
  std::cout << "\nDifference between coefficients computed by the two methods\n"
            << (gamma1 - gamma2) << std::endl;

  std::cout << "\nL2-Norms: " << std::sqrt(err1.sum()) << " "
            << std::sqrt(err2.sum()) << std::endl;

  std::cout << "\nDifference of the methods in L2-Norms: "
            << (std::sqrt(err1.sum()) - std::sqrt(err2.sum())) << std::endl;
}
