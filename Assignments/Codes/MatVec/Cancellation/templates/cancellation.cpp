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
  // TODO: Compute approximation of the derivative of sin(x)

  // Plot
  mgl::Figure fig;
  fig.setlog(true, true);
  fig.legend();
  fig.title("Error of approximation of f'(x_0)");
  fig.xlabel("h");
  fig.ylabel("| f'(x_0) - g_i(x_0, h) |");
  // TODO: Plot the error
  fig.save("error_cancellation.eps");
}
/* SAM_LISTING_END_1 */
