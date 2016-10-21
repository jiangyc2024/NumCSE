# include <iostream>
# include "fftreal.hpp"

int main() {
  Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(10, 0, 1);
  Eigen::VectorXcd c;

  fftreal(y, c);
  std::cout << "c:\n" << c << "\n";
  return 0;
}
