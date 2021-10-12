#include "naturalcubicspline.hpp"
#include "notaknotcubicspline.hpp"
#include <iostream>

int main(int /*argc*/, char ** /*argv*/) {
  constexpr double a = 0.0;
  constexpr double b = 4.0;
  constexpr int N = 5;
  Eigen::ArrayXd t{Eigen::ArrayXd::LinSpaced(N, a, b).pow(2)};
  Eigen::ArrayXd y{t.cos()};
  NaturalCubicSpline spline(t.matrix(), y.matrix());

  {
    const long n = t.size() - 1;
    constexpr double delta = 1.0E-8;
    std::cout << "Test of evalutation of cubic spline" << std::endl;
    std::cout << "t[0] = " << t[0] << ": y = " << y[0]
              << " <-> s(tp) = " << spline.eval(t[0] + delta) << std::endl;
    for (int i = 1; i < n; ++i) {
      const double val = y[i];
      const double knot = t[i];
      // Evaluating spline from left and right of the knots at a distance delta
      const double tm = knot - delta;
      const double tp = knot + delta;
      std::cout << "t[" << i << "] = " << knot
                << ": s(tm) = " << spline.eval(tm) << " <-> y = " << val
                << " <-> s(tp) = " << spline.eval(tp) << std::endl;
    }
    std::cout << "t[n] = " << t[n] << ": s(tm) = " << spline.eval(t[n] - delta)
              << " <-> y = " << y[n] << std::endl << std::endl;
  }

  // The below should be identical to what is printed for NaturalCubicSpline
  // above, as both splines match the data values at the knots. This will
  // already be the case before you start implementing NotAKnotCubicSpline,
  // since it is just a copy of the code for NaturalCubicSpline
  NotAKnotCubicSpline naks(t.matrix(), y.matrix());
  {
    const long n = t.size() - 1;
    constexpr double delta = 1.0E-8;
    std::cout << "Test of evalutation of NAK splines" << std::endl;
    std::cout << "t[0] = " << t[0] << ": y = " << y[0]
              << " <-> s(tp) = " << naks.eval(t[0] + delta) << std::endl;
    for (int i = 1; i < n; ++i) {
      const double val = y[i];
      const double knot = t[i];
      // Evaluating spline from left and right of the knots at a distance delta
      const double tm = knot - delta;
      const double tp = knot + delta;
      std::cout << "t[" << i << "] = " << knot << ": s(tm) = " << naks.eval(tm)
                << " <-> y = " << val << " <-> s(tp) = " << naks.eval(tp)
                << std::endl;
    }
    std::cout << "t[n] = " << t[n] << ": s(tm) = " << naks.eval(t[n] - delta)
              << " <-> y = " << y[n] << std::endl << std::endl;
  }

  {
    std::cout
        << "Test: A NAK spline should interpolate a cubic polynomial exactly"
        << std::endl;
    Eigen::ArrayXd p{t.pow(3)};
    NotAKnotCubicSpline nakp(t.matrix(), p.matrix());
    Eigen::ArrayXd x{Eigen::ArrayXd::LinSpaced(100, a + 1.0E-8, b - 1.0E-8)};
    // Norm of difference in spline and cubic polynomial evaluated for pts. in x
    const double diff =
        (x.unaryExpr([&nakp](double pt) { return nakp.eval(pt); }) - x.pow(3))
            .matrix()
            .norm();
    std::cout << "Norm of mismatch for cubic polynomial= " << diff << std::endl;
  }
  return 0;
}
