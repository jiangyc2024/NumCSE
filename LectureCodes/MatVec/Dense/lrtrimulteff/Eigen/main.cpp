#define NDEBUG true
#include <Eigen/Dense>
#include <figure/figure.hpp>
#include <iostream>
#include <numeric>
#include <cmath>
#include <cassert>
using namespace Eigen;
#include "lrtrimulteff.hpp"
int main(){
  int n = 3;
  MatrixXd A(n,n);
  A << 1,2,3,4,5,6,7,8,9;
  MatrixXd B(n,n);
  B << 9,8,7,6,5,4,3,2,1;
  VectorXd x(n);
  x << 4,5,6;
  VectorXd y(n);
  lrtrimulteff(A, B, x, y);
  std::cout << y << std::endl;
  return 0;
}
