# include <Eigen/Dense>
# include "figure.hpp"

int main()
{
  const unsigned int N = 1000;
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N, -5, 5);
  Eigen::VectorXd y = (1./(1 + x.array()*x.array())).matrix();

  Figure fig;
  fig.plot(x, y, "b");
  fig.save("runge.eps");

  return 0;
}
