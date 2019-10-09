//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <iostream>
#include <cmath>

#include <Eigen/Dense>
#include "linfit.hpp"

#include "matplotlibcpp.h"

using namespace Eigen;
namespace plt = matplotlibcpp;


int main() {
  VectorXd t = VectorXd::LinSpaced(10, 0.1, 1.0);
  MatrixXd A = make_A(t);
  
  std::cout << "\nThe matrix A is defined by \n" << A << std::endl;
  
  VectorXd f(10);
  f << 100. , 34. , 17. , 12. , 9. , 6. , 5. , 4. , 4. , 2.;
  
  
  VectorXd gamma;
  VectorXd tl = VectorXd::LinSpaced(91, 0.1, 1.0);
  VectorXd yl;
  VectorXd err;
  
  
  // approximate function coefficients using either method
  std::cout << "Enter \"1\" to only test data_fit_normal().\n"
              << "Enter \"2\" to only test data_fit_qr().\n";
  int ans=0;
  std::cin >> ans;
  
  switch(ans) {
    case 1: gamma = data_fit_normal(A, f);
            std::cout<< "Using data_fit_normal: gamma =\n"
                        << gamma << std::endl;
            fit_plot( gamma, tl, yl );
            error_plot( gamma, A, f, err );
            
            // plot data points and fitted function
            plt::figure();
            plt::semilogy(tl, yl, {{"label", "normal equation"}} );
            plt::semilogy(t, f, "o",{{"label","data"}} );
            plt::title("Fitted function");
            plt::xlabel("t");
            plt::ylabel("y");
            plt::legend("best");
            plt::savefig("./cx_out/fit.png");
            // plot fitting errors
            plt::figure();
            plt::semilogy(t, err, "o" , {{"label", "normal equation"}} );
            plt::title("Fit error");
            plt::xlabel("t");
            plt::ylabel("|y-f(t)|^2");
            plt::legend("best");
            plt::savefig("./cx_out/error.png");
            break;
    case 2: gamma = data_fit_qr(A, f);
            std::cout<< "Using data_fit_qr: gamma =\n"
                        << gamma << std::endl;
            fit_plot( gamma, tl, yl );
            error_plot( gamma, A, f, err );
            
            // plot data points and fitted function
            plt::figure();
            plt::semilogy(tl, yl, {{"label", "qr fitting"}} );
            plt::semilogy(t, f, "o",{{"label","data"}} );
            plt::title("Fitted function");
            plt::xlabel("t");
            plt::ylabel("y");
            plt::legend("best");
            plt::savefig("./cx_out/fit.png");
            // plot fitting errors
            plt::figure();
            plt::semilogy(t, err, "o" , {{"label", "qr fitting"}} );
            plt::title("Fit error");
            plt::xlabel("t");
            plt::ylabel("|y-f(t)|^2");
            plt::legend("best");
            plt::savefig("./cx_out/error.png");
            break;
    default: return 0;
  }
}
