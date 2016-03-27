# include "./intpolyval_lag.hpp"
# include <figure/figure.hpp>

int main () {
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(6, -1, 1),
                  y = ((5*t).array().sin()*(-t).array().exp()).matrix(),
                  x = Eigen::VectorXd::LinSpaced(200, -1, 1),
                  p;

  intpolyval_lag(t, y, x, p);
  mgl::Figure fig;
  fig.title("Intp. with Lagrange polynomials");
  fig.plot(t, y, " r*").label("Data");
  fig.plot(x, p).label("Interpolant");
  fig.fplot("sin(5*x)*exp(-x)").label("f(x)");
  fig.legend();

  fig.save("intp");

  return 0;
}
