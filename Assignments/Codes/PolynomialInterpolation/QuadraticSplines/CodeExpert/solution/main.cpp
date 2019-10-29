//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <iostream>

#include "quadsplines.hpp"

using namespace Eigen;


int main(void) {
  VectorXd t(3);
  t << .1, .2, .5 ;
  std::pair<VectorXd, VectorXd> p = increments(t);
  std::cout << p.first << std::endl << p.second << std::endl;
  
  VectorXd y(4);
  y << 1, 2, 1, 2;
  
  VectorXd c = compute_c(t, y);
  std::cout << "c = " << c << std::endl << std::endl; 
  
  Vector3d x(0,.1,.9) ;
  
  VectorXd fval = quadspline(t, y, x);
  std::cout << "fval = " << fval << std::endl << std::endl; 
  plotquadspline("splines");
  
  
  std::vector<double> Err = qsp_error(8);
  std::cout << "Error values: ";
  for (unsigned int i = 0; i < Err.size(); ++i) {
    std::cout << Err[i] << "    " ;
  }
  
  return 0;
}
