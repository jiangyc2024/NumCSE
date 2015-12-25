// this is a short example code showing how the library can be used
# include <iostream>
# include "figure.hpp"

int main()
{
  Eigen::VectorXd x(5), y(5);
  x << 1,2,3,4,5;
  y << 1,2,1,2,5;
  std::vector<double> v = {3, 0.5, 6}; // note: must be vector of doubles!
  std::vector<double> w = {0.1, 1, 3};

  Figure fig;
  fig.title("Testplot");
  fig.setlog(true, true);
  fig.plot(x, y, "r#+", "My data");
  fig.plot(v, w, "b*", "Other data");
  fig.grid();
  fig.xlabel("x-Axis");
  fig.legend();
  fig.save("dataplot.eps");

  return 0;
}
