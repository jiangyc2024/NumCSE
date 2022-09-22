///
/// Minimal runner for (0-3)
///

#include <iostream>

#include "MatrixReduce.hpp"

int main() {
  std::cout << "\nTest of average():\n";
  Eigen::Matrix3d A = Eigen::Matrix3d::Identity();
  std::cout << average(A) << "\n\n";
  std::cout << "Test of percent_zero():\n";
  std::cout << percent_zero(A) << "\n\n";

  std::cout << "Test of has_zero_column():\n";
  Eigen::MatrixXd B = Eigen::MatrixXd::Identity(4, 5);
  std::cout << has_zero_column(B) << " " << has_zero_column(B.transpose())
            << "\n\n";

  std::cout << "Test of columns_sum_to_zero():\n";
  std::srand(5);  // So that random behaviour is predictable.
  Eigen::Matrix3d C = Eigen::Matrix3d::Random() + Eigen::Matrix3d::Constant(1);
  std::cout << columns_sum_to_zero(C) << "\n\n";

  return 0;
}
