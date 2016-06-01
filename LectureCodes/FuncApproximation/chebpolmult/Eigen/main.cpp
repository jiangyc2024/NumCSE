# include <iostream>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "chebpolmult.hpp"

int main() {
 
  Eigen::RowVectorXd x = Eigen::RowVectorXd::LinSpaced(100, -1, 1);
  Eigen::MatrixXd V;
  unsigned d = 5;
  chebpolmult(d, x, V);

  mgl::Figure fig;
  for (unsigned r = 0; r < d + 1; ++r) {
    fig.plot(x, V.row(r)).label("k = " + std::to_string(r));
  }
  fig.title("First " + std::to_string(d) + " cheb. polynomials");
  fig.xlabel("x");
  fig.ylabel("T_k(x)");
  fig.ranges(-1,1,-1.1,1.1);
  fig.legend();
  fig.save("chebs.eps");
  

  return 0;
}
