
#include <iostream>
#include "LV.hpp"

using namespace Eigen;


int main(){
  
  // The test uses the input u0=2.8, v0=1.5, T=2
  std::pair<Vector2d, Matrix2d> PaW = PhiAndW(2.8, 1.5, 2);
  std::cout << "Test of PhiAndW():\nPhi = "
            << PaW.first.transpose()
            << "\nW = \n"
            << PaW.second << "\n\n";
  
  Vector2d y = findInitCond();
  std::cout << "Test of findInitCond():\ny = "
            << y.transpose() << "\n";
  
  return 0;
}