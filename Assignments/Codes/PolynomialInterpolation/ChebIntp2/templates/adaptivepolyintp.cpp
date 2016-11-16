# include <cmath>
# include <figure/figure.hpp>

#include "intpolyval.hpp"


/*!
 * \brief adaptivepolyintp
 * \tparam Function
 * \param f
 * \param a
 * \param b
 * \param tol
 * \param N
 * \param tRes
 * \param errRes
 */
/* SAM_LISTING_BEGIN_1 */
template <class Function>
void adaptivepolyintp(const Function& f, const double a, const double b,
                      const double tol, const unsigned N,
                      Eigen::VectorXd& adaptive_nodes,
                      /* Ignore the following line for part a) */
                      Eigen::VectorXd& error_vs_step_no) {
    // Generate sampling points and evaluate $f$ there
    Eigen::VectorXd sampling_points =
            Eigen::VectorXd::LinSpaced(N, a, b),
                    fvals_at_sampling_points =
            sampling_points.unaryExpr(f);

    // TODO: perform adaptive interpolation and return nodes
    // (and error for part b)
}
/* SAM_LISTING_END_1 */


int main() {
  /* SAM_LISTING_BEGIN_2 */
  // Declare test functions
  auto f1 = [](double t) { return std::sin(std::exp(2*t)); };
  auto f2 = [](double t) { return std::sqrt(t)/(1 + 16*t*t); };

  // Test interval
  const double a = 0, b = 1;

  // Get interpolation nodes and print runtimes
  const unsigned N = 1000; // no. of sampling points
  const double tol = 1e-6; // tolerance
  Eigen::VectorXd tf1, tf2, // nodes for f1 resp. f2
                  ef1, ef2; // errors for f1 resp. f2

  // TODO: do adaptive interpolation for f1 and f2

  // Plot
  mgl::Figure fig;
  fig.title("Error VS step");
  fig.setlog(false, true);
  fig.xlabel("No. of interpolation nodes");
  fig.ylabel("max |f(t) - I_Tf(t)|");
  fig.plot(ef1, "r").label("f_1(t) = sin(e^{2t})");
  fig.plot(ef2, "b").label("f_2(t) = \\sqrt{t}/(1 + 16t^2)");
  fig.legend();
  fig.save("cvg");
  /* SAM_LISTING_END_2 */

  return 0;
}
