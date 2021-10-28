#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

#include "multAmin.hpp"

int main() {
  constexpr unsigned int M = 10;
  Eigen::VectorXd xa = Eigen::VectorXd::Random(M);
  Eigen::VectorXd ys, yf;

  std::cout << "\nEnter \"0\" to test all functions.\n"
            << "Enter \"1\" to only test multAmin().\n"
            << "Enter \"2\" to only test multAmin_runtime().\n"
            << "Enter \"3\" to only test multABunitv().\n";

  int ans = 0;
  std::cin >> ans;
  switch (ans) {
    case 0:  // Testing correctness of the code
      multAmin(xa, yf);
      multAminSlow(xa, ys);

      // Checking error. Error should be small!
      std::cout << "\nChecking error and comparing runtimes of multAmin()"
                << " and multAminSlow().\nThe error should be small!\n";

      std::cout << "\nThe error is: ||y_slow-y_fast|| = " << (ys - yf).norm()
                << std::endl
                << "\n";

      multAmin_runtime();
      std::cout << std::defaultfloat;

      // Testing multABunitv()
      std::cout << "\nTesting multABunitv()\n";
      multABunitv();
      break;
    case 1:  // Testing correctness of the code
      multAmin(xa, yf);
      multAminSlow(xa, ys);

      std::cout << "\nError between result of multAmin() and multAminSlow().\n"
                << "The error should be small.\n";

      // Error should be small
      std::cout << "\n||y_slow-y_fast|| = " << (ys - yf).norm() << std::endl;
      break;
    case 2:
      multAmin_runtime();
      break;
    case 3:
      Eigen::MatrixXd C = multABunitv();
      break;
  }
}
