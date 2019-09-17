
#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>

#include "multAmin.hpp"

using namespace Eigen;


int main(void) {
  unsigned int M = 10;
  VectorXd xa = VectorXd::Random(M);
  VectorXd ys, yf;
  
  std::cout << "\nEnter \"0\" to test all functions.\n"
              << "Enter \"1\" to only test multAmin().\n"
              << "Enter \"2\" to only test multAmin_runtime().\n"
              << "Enter \"3\" to only test multABunitv().\n";
  
  int ans=0;
  std::cin >> ans;
  switch(ans){
    case 0: // Testing correctness of the code
            multAmin(xa, yf);
            multAminSlow(xa, ys);
            // Error should be small
            std::cout << "||y_slow-y_fast|| = " << (ys - yf).norm() << std::endl;
            multAmin_runtime();
            std::cout << std::defaultfloat;
            multABunitv();
            break;
    case 1: // Testing correctness of the code
            multAmin(xa, yf);
            multAminSlow(xa, ys);
            // Error should be small
            std::cout << "||y_slow-y_fast|| = " << (ys - yf).norm() << std::endl;
            break;
    case 2: multAmin_runtime();
            break;
    case 3: MatrixXd C = multABunitv();
            break;
  }    
    
}
