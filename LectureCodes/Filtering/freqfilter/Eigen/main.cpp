#include "freqfilter.hpp"
#include <iostream>

using Eigen::VectorXd;

int main() { //NOLINT(bugprone-exception-escape)
  const int n = 10;
  const int k = 2;
  const VectorXd y = Eigen::ArrayXd::LinSpaced(n, 0, M_PI).cos().matrix();
  VectorXd l;
  VectorXd h;
  freqfilter::freqfilter(y, k, l, h);
  std::cout << "Low:\n" << l.transpose() << "\n"
            << "High:\n" << h.transpose() << "\n";
 return 0;
}
