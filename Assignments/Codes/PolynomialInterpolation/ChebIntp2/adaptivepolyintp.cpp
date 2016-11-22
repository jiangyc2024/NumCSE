#include <cmath>
#include <mgl2/mgl.h>
#include <figure/figure.hpp>

#include "intpolyval.hpp"

#if INTERNAL
#include <timer.h>
#endif // INTERNAL

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

# if SOLUTION

    // Approximate $\max |f(x)|$
    double fmax = fvals_at_sampling_points.cwiseAbs().maxCoeff();

    // Adaptive mesh (initial node)
    std::vector<double> t { (a+b)/2. }, // Set of interpolation nodes
                        y { static_cast<double>(f((a+b)/2.)) }, // Values at nodes
                        errors;         // Error at nodes


    for (int i = 0; i < N; ++i) {
        // *** Step 1: interpolate with current nodes
        //   need to convert std::vector to
        //   Eigen::VectorXd to use the function interpoyval
        Eigen::Map<Eigen::VectorXd> te(t.data(), t.size());
        Eigen::Map<Eigen::VectorXd> ye(y.data(), y.size());
        Eigen::VectorXd intpolyvals_at_sampling_points;
        intpolyval(te, ye, sampling_points, intpolyvals_at_sampling_points);

        // *** Step 2: find node where error is the largest
        Eigen::VectorXd err =
                (fvals_at_sampling_points - intpolyvals_at_sampling_points)
                .cwiseAbs();
        double max = 0; int idx = 0;

        // We use an Eigen "Visitor"
        /*
        https://eigen.tuxfamily.org/dox/group__TutorialReductionsVisitorsBroadcasting.html
        */
        max = err.maxCoeff(&idx);

        // Alternative:
        /*
        for (int j = 0; j < err.size(); ++j) {
            if (err(j) > max) {
              max = err(j); idx = j;
            }
        }
        */

        // Step 3: check termination criteria
        if (max < tol * fmax) {
            // If terminate save results in correct variables
            adaptive_nodes = te;
            Eigen::Map<Eigen::VectorXd> tmp(errors.data(), errors.size());
            error_vs_step_no = tmp;
            return;
        }

        // Step 4 (part b): add this node to our set of nodes and save error
        errors.push_back(max);
        t.push_back(sampling_points(idx));
        y.push_back(fvals_at_sampling_points(idx));
    }
    std::cerr << "Desired accuracy could not be reached."
              << std::endl;
    adaptive_nodes = sampling_points; // return all sampling points

#else // TEMPLATE
    // TODO: perform adaptive interpolation and return nodes
    // (and error for part b)
#endif // TEMPLATE
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

#if SOLUTION
  adaptivepolyintp(f1, a, b, tol, N, tf1, ef1);
  adaptivepolyintp(f2, a, b, tol, N, tf2, ef2);
#else // TEMPLATE
  // TODO: do adaptive interpolation for f1 and f2
#endif // TEMPLATE

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