#ifndef GAUSSCONVERGENCE_HPP
#define GAUSSCONVERGENCE_HPP

#include <Eigen/Dense>
#include <cmath>
#include <vector>

#include "gaussquad.hpp"

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

/*!
 * \brief integrate Compute integral given quadrature rule.
 * Compute the integral $\int_a^b f(x) dx \approx \sum_{i=1}^n w_i f(c_i)$
 * for given quadrature rule $\{(w_i, x_i)\}_{i=1}^n$
 * \tparam Function A function object with operator()(double)
 * \param qr A QuadRule object passing nodes and weights
 * \param f A function handle passing the integrand
 * \return Value of integral of $f$ using quadrature rule Q
 */
/* SAM_LISTING_BEGIN_0 */
template <class Function>
double integrate(const QuadRule& qr, const Function& f) {
  double I = 0;
  // Compute weighted sum of function values according to
  // the concept of a generic quadrature formula 
  for (unsigned i = 0; i < qr.weights_.size(); ++i) {
    I += qr.weights_(i) * f(qr.nodes_(i));
  }
  return I;
}
/* SAM_LISTING_END_0 */

/*!
 * \brief gaussConv Approximte the integral $\int_{-1}^1 \arcsin(x) f(x) dx$
 *
 * \tparam Function A function object with operator()(double)
 * \param fh Will pass the integrand
 * \param I_ex Exact value of integral
 * \return Value of integral
 */
/* SAM_LISTING_BEGIN_1 */
template <class Function>
double gaussConv(const Function& fh, const double I_ex, const unsigned N) {
  std::vector<double> evals,  // used to save no.\ of quad nodes
      error,                  // used to save the error
      evals_ref;  // used to plot reference line of order of convergence

  double I_N = 0.;  // best approx. of integral
  plt::figure();

  // TO DO (8-1.a): write code that plots the no. of nodes versus the error
  //                using matplotlibcpp and returns the best approximation
  // START

  // Build integrand
  auto f = [fh](double x) { return std::asin(x) * fh(x); };
  
  double I;
  for (unsigned n = 1; n <= N; ++n) {
    QuadRule qr = gaussquad(n);  // Create quadrature rule

    I = integrate(qr, f);  // Compute integral

    evals.push_back(n);  // Save no. of quadrature nodes
    evals_ref.push_back(
        std::pow(n, -3));  // used for plotting to show convergence rate
    error.push_back(std::abs(I - I_ex));  // Save error
  }
  // save best approx. of integral
  I_N = I;

  // Create convergence plots
  plt::title("Gauss quadrature convergence");
  plt::loglog(evals, error, " +r", {{"label", "Error"}});  // plot error
  plt::loglog(evals, evals_ref, "k--",
              {{"label", "O(n$^{-3}$)"}});  // reference line
  plt::xlabel("No. of quadrature nodes");
  plt::ylabel("|Error|");
  plt::legend("best");
  // END
  
  plt::savefig("./cx_out/GaussConv.png");
  // return the best approximation
  return I_N;
}
/* SAM_LISTING_END_1 */

/*!
 * \brief gaussConv Approximte the integral $\int_{-1}^1 \arcsin(x) f(x) dx$
 * Ensures that convergenge is expoenential using appropriate transformation,
 * provided $f$ is a smooth function.
 * \tparam Function A function object with operator()(double)
 * \param fh Will pass the integrand
 * \param I_ex Exact value of integral
 * \return Value of integral
 */
/* SAM_LISTING_BEGIN_2 */
template <class Function>
double gaussConvCV(const Function& f, const double I_ex, const unsigned N) {
  std::vector<double> evals,  // Used to save no. of quad nodes
      error;                  // Used to save the error

  double I_N = 0.;  // best approx. of integral
  plt::figure();

  // TO DO (8-1.d): complete to perform Gauss-Legendre quadrature with the
  // variable transform of the integral derived in 8-1.c)
  //                write code that plots the no. of nodes versus the error
  //                using matplotlibcpp and returns the best approximation
  // START

  // Lambda function providing \cor{transformed integrand}
  auto g = [f](double x) { return x * f(std::sin(x)) * std::cos(x); };

  double I;
  for (unsigned n = 1; n <= N; ++n) {
    QuadRule qr = gaussquad(n);  // Obtain quadrature rule

    // Transform nodes and weights to new interval
    Eigen::VectorXd w = qr.weights_ * M_PI / 2;
    Eigen::VectorXd c =
        (-M_PI / 2 + M_PI / 2 * (qr.nodes_.array() + 1)).matrix();

    // Evaluate $g$ at quadrature nodes $c$
    Eigen::VectorXd gc = c.unaryExpr(g);

    // Same as $I = \sum_{i=1}^n w_i g(c_i)$
    I = w.dot(gc);

    evals.push_back(n);                 // Save no. of quadrature nodes
    error.push_back(std::abs(I - I_ex));  // Save error
  }
  
  // save approximation of integral computed with highest-order rule
  I_N = I;

  // Create convergence plots
  plt::title("Gauss quadrature convergence");
  plt::semilogy(evals, error, " +r", {{"label", "Error"}});  // plot error
  plt::xlabel("No. of quadrature nodes");
  plt::ylabel("|Error|");
  plt::legend("best");
  // END

  plt::savefig("./cx_out/GaussConvCV.png");
  // return the best approximation
  return I_N;
}
/* SAM_LISTING_END_2 */

#endif

