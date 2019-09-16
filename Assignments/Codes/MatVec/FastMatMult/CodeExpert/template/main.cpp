#include <Eigen/Dense>
#include <iostream>

#include "strassen.hpp"

using namespace Eigen;

int main() {
  
  double tolerance = 1e-9;
  int n = 128;
  std::srand(5);
  MatrixXd A = MatrixXd::Random(n,n);
  std::cout << "Works for identity matrix: "
            << ((strassenMatMult(A, MatrixXd::Identity(n,n)) - A).norm()<tolerance)
            << " "
            << ((strassenMatMult(MatrixXd::Identity(n,n),A) - A).norm()<tolerance)
            << std::endl;
  std::cout << "Works for two random matrices: "
            <<(test_strassen() < tolerance) << std::endl;
  
  std::cout << "Output of test_strassen()" << test_strassen() << std::endl;
  time_strassen();
}
