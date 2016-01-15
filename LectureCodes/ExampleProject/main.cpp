# include <Eigen/Dense>
# include "figure.hpp"

int main()
{
  const unsigned int N = 1000;
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N, -5, 5);
  Eigen::VectorXd y = (1./(1 + x.array()*x.array())).matrix();

  Figure fig;
  fig.plot(x, y, "b", "1/\\(1 + x^2)^{-1}"); // Latex typesetting is supported when typing '\\' before the commands
  fig.legend();
  fig.save("runge.eps");

  return 0;
}
