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
  
  VectorXd gamma1 = data_fit_normal(t, f);
  VectorXd gamma2 = data_fit_qr(t, f);
  
  VectorXd tl = VectorXd::LinSpaced(91, 0.1, 1.0);
  VectorXd yl1, yl2;
  VectorXd err1, err2;
  
  error_plot( gamma1, A, f, err1 );
  error_plot( gamma2, A, f, err2 );
  // approximate function coefficients using both methods
  std::cout << "Enter \"1\" to plot data_fit_normal().\n"
              << "Enter \"2\" to plot data_fit_qr().\n"
              << "Enter another key to exit.\n";
  int ans=0;
  std::cin >> ans;
  
  switch(ans) {
    case 1: 
            std::cout<< "Using data_fit_normal: gamma =\n"
                        << gamma1 << std::endl;
            fit_plot( gamma1, tl, yl1 );
            
            // plot data points and fitted function
            plt::figure();
            plt::semilogy(tl, yl1, {{"label", "normal equation"}} );
            plt::semilogy(t, f, "o",{{"label","data"}} );
            plt::title("Fitted function");
            plt::xlabel("t");
            plt::ylabel("y");
            plt::legend("best");
            plt::savefig("./cx_out/fit1.png");
            // plot fitting errors
            plt::figure();
            plt::semilogy(t, err1, "o" , {{"label", "normal equation"}} );
            plt::title("Fit error");
            plt::xlabel("t");
            plt::ylabel("|y-f(t)|^2");
            plt::legend("best");
            plt::savefig("./cx_out/error1.png");
            break;
    case 2: 
            std::cout<< "Using data_fit_qr: gamma =\n"
                        << gamma2 << std::endl;
            fit_plot( gamma2, tl, yl2 );
            
            // plot data points and fitted function
            plt::figure();
            plt::semilogy(tl, yl2, {{"label", "qr fitting"}} );
            plt::semilogy(t, f, "o",{{"label","data"}} );
            plt::title("Fitted function");
            plt::xlabel("t");
            plt::ylabel("y");
            plt::legend("best");
            plt::savefig("./cx_out/fit2.png");
            // plot fitting errors
            plt::figure();
            plt::semilogy(t, err2, "o" , {{"label", "qr fitting"}} );
            plt::title("Fit error");
            plt::xlabel("t");
            plt::ylabel("|y-f(t)|^2");
            plt::legend("best");
            plt::savefig("./cx_out/error2.png");
            break;
    default: return 0;
  }
  // show that the different methods don't give the same result
  std::cout << "\nDifference between two methods\n" 
              << (gamma1 - gamma2) << std::endl;
  std::cout << "\nL2-Norms: " << std::sqrt(err1.sum()) 
              << " " << std::sqrt(err2.sum()) << std::endl;
  std::cout << "Difference of the methods in L2-Norms: " 
              << (std::sqrt(err1.sum()) - std::sqrt(err2.sum())) << std::endl;
}
