# include "./freqfilter.hpp"
# include <iostream>

int main() {
  const int n = 10,
            k = 2;
  VectorXd y = Eigen::ArrayXd::LinSpaced(n, 0, M_PI).cos().matrix();
  VectorXd l, h;
  freqfilter(y, k, l, h);
  std::cout << "Low:\n" << l.transpose() << "\n"
            << "High:\n" << h.transpose() << "\n";
 return 0;
}
