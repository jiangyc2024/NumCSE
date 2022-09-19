///
/// Minimal runner for (0-2)
///

#include <iostream>

#include "MatrixBlocks.hpp"

int main() {
  std::cout << "\nTest of zero_row_col():\n";
  Eigen::Matrix3d A = Eigen::Matrix3d::Constant(-1);
  std::cout << zero_row_col(A, 0, 1) << "\n\n";

  std::cout << "Test of swap_left_right_blocks():\n";
  Eigen::MatrixXd B = Eigen::MatrixXd::Identity(4, 3);
  std::cout << swap_left_right_blocks(B, 2) << "\n\n";

  std::cout << "Test of tridiagonal():\n";
  std::cout << tridiagonal(4, -1, 2, -1) << "\n\n";
  return 0;
}
