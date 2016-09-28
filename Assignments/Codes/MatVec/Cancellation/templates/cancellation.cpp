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
