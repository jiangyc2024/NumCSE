# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "feval.hpp"
# include "polyfit.hpp"
# include "polyval.hpp"


int main() {
#pragma
  auto f = [](double x) { return 1./(1 + x*x); };
  // sampling points for approximate maximum norm
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(1000, -5, 5), 
                 fx = feval(f, x); // evaluate f at x

  std::vector<double> err;
  for (int d = 1; d <= 20; ++d) {
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(d + 1, -5, 5),
                   ft = feval(f, t);

    Eigen::VectorXd p = polyfit(t, ft, d),
                    y = polyval(p, x);
    err.push_back( (y - fx).cwiseAbs().maxCoeff() );
  }
#pragma

  mgl::Figure fig;
  fig.setlog(false, true);
  fig.plot(err, "r--+");

  fig.xlabel("degree d");
  fig.ylabel("interpolation error (maximum norm)");
  fig.save("rungeerrmax.eps");

  return 0;
}
