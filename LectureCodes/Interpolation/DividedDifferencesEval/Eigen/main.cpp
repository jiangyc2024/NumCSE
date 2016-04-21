# include "./evaldivdiff.hpp"
# include <figure/figure.hpp>
# include <sstream>
# include <string>

int main() {
  /*
   * testing divdiff code for (0) runge function and (1) sine 
   */

  mgl::Figure f0, f1;
  f0.title("f(t) = \\frac{1}{1 + t^2}");
  f0.ranges(-5, 5, -0.2, 1.2);
  f0.legend();
  f0.fplot("1/(1 + x^2)").label("Exact");

  f1.title("f(t) = sin(x)");
  f1.ranges(0, 2*M_PI, -1.2, 1.2);
  f1.legend();
  f1.fplot("sin(x)").label("Exact");

  for (unsigned n = 4; n <= 10; n += 2) {
    // (0) : Runge function 
    VectorXd t0 = VectorXd::LinSpaced(n, -5, 5),
             y0 = ( 1./(1 + t0.array()*t0.array()) ).matrix();

    // (1) : sine
    VectorXd t1 = VectorXd::LinSpaced(n, 0, 2*M_PI),
             y1 = t1.array().sin().matrix();

    // evaluating
    VectorXd x0 = VectorXd::LinSpaced(100, -5, 5),
             x1 = VectorXd::LinSpaced(100, 0, 2*M_PI);
    VectorXd p0, p1;
    evaldivdiff(t0, y0, x0, p0);
    evaldivdiff(t1, y1, x1, p1);

    // plot
    std::stringstream ss;
    ss << n;
    f0.plot(x0, p0).label("n = " + ss.str());
    f1.plot(x1, p1).label("n = " + ss.str());
  }

  // save plots
  f0.save("dd_runge");
  f1.save("dd_sine");

  return 0;
}
