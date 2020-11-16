#ifndef ADAPTINTERP_HPP
#define ADAPTINTERP_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>
#include "intpolyval.hpp"
#include "matplotlibcpp.h"

using namespace Eigen;
namespace plt = matplotlibcpp;
/*!
 * \brief adaptivepolyintp
 * \tparam Function
 * \paramin f
 * \paramin a
 * \paramin b
 * \paramin tol
 * \paramin N
 * \param errortab, pointer to vector
 */
/* SAM_LISTING_BEGIN_1 */
template <class Function>
VectorXd adaptivepolyintp(Function&& f, double a, double b, double tol, int N,
                          std::vector<double>* errortab = nullptr) {
  // TO DO (7-3.a) : implement the greedy algorithm for adaptive interpolation.
  // Ignore the errortab part of this function for now.

  // Generate sampling points and evaluate $f$ there
  VectorXd sampling_points = VectorXd::LinSpaced(N, a, b),
           fvals_at_sampling_points = sampling_points.unaryExpr(f);
  // START
  // Approximate $\max |f(x)|$
  double maxf = fvals_at_sampling_points.cwiseAbs().maxCoeff();

  // Adaptive mesh (initial node)
  std::vector<double> t{(a + b) / 2.},          // Set of interpolation nodes
      y{static_cast<double>(f((a + b) / 2.))},  // Values at nodes
      errors;                                   // Error at nodes

  for (int i = 0; i < N; ++i) {
    // *** Step 1: interpolate with current nodes
    //   need to convert std::vector to
    //   Eigen::VectorXd to use the function intpolyval
    Map<VectorXd> te(t.data(), t.size());
    Map<VectorXd> ye(y.data(), y.size());
    VectorXd intpolyvals_at_sampling_points =
        intpolyval(te, ye, sampling_points);

    // *** Step 2: find node where error is the largest
    VectorXd err =
        (fvals_at_sampling_points - intpolyvals_at_sampling_points).cwiseAbs();
    // We use an Eigen "Visitor"
    // https://eigen.tuxfamily.org/dox/group__TutorialReductionsVisitorsBroadcasting.html

    double max = 0;
    int idx = 0;
    max = err.maxCoeff(&idx);

    // Step 3: check termination criteria
    if (max < tol * maxf) {
      return te;
    }
    // Step 4: add this node to our set of nodes
    t.push_back(sampling_points(idx));
    y.push_back(fvals_at_sampling_points(idx));
    // TO DO (7-3.b): save error
    if (errortab != nullptr) {
      errortab->push_back(max);
    }
  }
  std::cerr << "Desired accuracy could not be reached." << std::endl;
  // END
  return sampling_points;  // return all sampling points
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void plotInterpolationError(void) {
  // Declare test functions
  auto f1 = [](double t) { return std::sin(std::exp(2 * t)); };
  auto f2 = [](double t) { return std::sqrt(t) / (1 + 16 * t * t); };
  // TO DO (7-3.c): generate the plots of error vs number of nodes for f1, f2

  // START
  // Test interval
  const double a = 0, b = 1;
  // Get interpolation nodes
  const unsigned N = 1000;  // no. of sampling points
  const double tol = 1e-6;  // tolerance

  VectorXd tf1, tf2;             // nodes for f1 resp. f2
  std::vector<double> ef1, ef2;  // errors for f1 resp. f2

  tf1 = adaptivepolyintp(f1, a, b, tol, N, &ef1);
  tf2 = adaptivepolyintp(f2, a, b, tol, N, &ef2);
  VectorXd n1 = VectorXd::LinSpaced(ef1.size(), 1, ef1.size());
  VectorXd n2 = VectorXd::LinSpaced(ef2.size(), 1, ef2.size());
  // Plot
  plt::figure();
  plt::title("Error VS step");
  plt::xlabel("No. of interpolation nodes");
  plt::ylabel("max |f(t) - I_Tf(t)|");
  plt::semilogy(n1, ef1, "ro", {{"label", "$f_1(t) = sin(e^{2t})$"}});
  plt::semilogy(n2, ef2, "bo", {{"label", "$f_2(t) = \\sqrt{t}/(1 + 16t^2)$"}});
  plt::legend();
  plt::savefig("./cx_out/intperrplot.png");
  // END
}
/* SAM_LISTING_END_2 */
#endif
