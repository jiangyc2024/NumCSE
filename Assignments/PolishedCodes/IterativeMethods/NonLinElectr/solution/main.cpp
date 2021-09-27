///
/// Minimal runner for (9-7)
///

#include <Eigen/Dense>
#include <iostream>

#include "nonlinear_circuit.hpp"

int main() {
  constexpr unsigned int n = 20;
  constexpr double alpha = 8;
  constexpr double beta = 1;
  
  Eigen::VectorXd Uin = Eigen::VectorXd::LinSpaced(n, 0, 20);
  Eigen::VectorXd Uout(n);
  circuit(alpha, beta, Uin, Uout);
  
  // Display the solution
  std::cout << "The solutions are:\n" << Uout << std::endl;
  
  // Display the differences of the solutions: the nonlinear effect can be seen from the fact that the vector is not constant
  Eigen::VectorXd diff(n - 1);
  for (unsigned int i = 0; i < n - 1; ++i) {
    diff(i) = Uout(i + 1) - Uout(i);
  }
  std::cout << "The differences are:\n" << diff << std::endl;
  
  plotU();
}
