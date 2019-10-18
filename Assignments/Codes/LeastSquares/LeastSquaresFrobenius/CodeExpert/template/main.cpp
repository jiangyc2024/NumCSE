//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <Eigen/Dense>
#include <iostream>

#include "matrixfit.hpp"

using namespace Eigen;

int main() {

  unsigned int n = 5;

  VectorXd z = VectorXd::LinSpaced(n,1,n),
           g(n) ;
  g << 5 , 4 , 3 , 2 , 1;
  MatrixXd Mstar = min_frob(z, g); // $M^*$
  
  std::cout << Mstar << std::endl;
  
  bool works = testMformula(100);
  if (works) {
    std::cout << "The formula gives the correct result"<< std::endl;
  }
  else {
    std::cout << "The two methods give different results"<< std::endl;
  }
 
}
