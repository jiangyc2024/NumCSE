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
  
  VectorXd d = compute_d(c,t);
  std::cout << "d = " << d << std::endl << std::endl; 
  
  Vector3d x(0,.4,.9) ;
  
  VectorXd fval = quadspline(t, y, x);
  std::cout << "fval = " << fval << std::endl << std::endl; 
  plotquadspline("splines");
  
  
  std::vector<double> Err = qsp_error(8);
  std::cout << "n \t" << "Error values: \n";
  std::cout << "------------------------ \n";
  for (unsigned int i = 0; i < Err.size(); ++i) {
    std::cout << std::pow(2,i+1) << "\t" << Err[i] << std::endl ;
  }
  
  return 0;
}
