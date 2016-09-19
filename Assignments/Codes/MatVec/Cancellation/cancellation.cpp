///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Till Ehrengruber <tille@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>

#include <figure/figure.hpp>

using namespace Eigen;

/* SAM_LISTING_BEGIN_1 */
int main() {
#if SOLUTION
  auto h = ArrayXd::LinSpaced(21, -20, 0).unaryExpr([] (const double i) {
    return std::pow(10, i);
  });
  auto x = ArrayXd::Constant(h.size(), 1.2);

  // Derivative
  ArrayXd g1 = (sin(x+h) - sin(x)) / h; // naive
  ArrayXd g2 = 2 * cos(x+0.5*h) * sin(0.5 * h) / h;
  ArrayXd ex = cos(x);
#else
  // TODO: Compute approximation of the derivative of sin(x)
#endif

  // Plot
  mgl::Figure fig;
  fig.setlog(true, true);
  fig.legend();
  fig.title("Error of approximation of f'(x_0)");
  fig.xlabel("h");
  fig.ylabel("| f'(x_0) - g_i(x_0, h) |");
#if SOLUTION
  fig.plot(h.matrix(), (g1-ex).abs().matrix()).label("g_1");
  fig.plot(h.matrix().tail(16), (g2-ex).abs().matrix().tail(16)).label("g_2");
  fig.plot(h.matrix(), h.matrix(), "h;").label("O(h)");
#else
  // TODO: Plot the error
#endif
  fig.save("error_cancellation.eps");
}
/* SAM_LISTING_END_1 */
