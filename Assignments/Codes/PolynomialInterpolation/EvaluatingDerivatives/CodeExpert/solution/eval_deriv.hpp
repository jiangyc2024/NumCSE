#include <complex>
#include <Eigen/Dense>
#include "timer.h"
#include "polyfit.hpp"
#include "polyval.hpp"
#include "matplotlibcpp.h"

using namespace Eigen;
namespace plt = matplotlibcpp;

/*
 * @brief Evaluate a polynomial and its derivative using Horner scheme
 * @param[in] vector c of size $n$, coefficients of the polynomial p
 * @param[in] double x, where the polynomial has to be evaluated
 * @param[out] pair containing p(x),p'(x)
 */
// TO DO (6-1.a): Implement the function evaldp().
// START
/* SAM_LISTING_BEGIN_0 */
template <typename CoeffVec>
std::pair<double, double> evaldp(const CoeffVec &c, const double x) {
  double px, dpx;
  int s = c.size(); // degree + 1
  // Standard Horner scheme
  px = c[0]; for (int i = 1; i < s; ++i) px = x * px + c[i];
  // Horner scheme for derivative  
  dpx = (s - 1) * c[0];
  for (int i = 1; i < s - 1; ++i) dpx = x * dpx + (s - i - 1) * c[i];

  return { px, dpx }; 
}
/* SAM_LISTING_END_0 */
// END


/*
 * @brief Evaluate a polynomial and its derivative using a naive implementation
 * @param[in] vector c of size $n$, coefficients of the polynomial p
 * @param[in] double x, where the polynomial has to be evaluated
 * @param[out] pair containing p(x),p'(x)
 */
// TO DO (6-1.b): Implement the function evaldp_naive().
// START
/* SAM_LISTING_BEGIN_1 */
template <typename CoeffVec>
std::pair<double, double> evaldp_naive(const CoeffVec &c, const double x) {
  std::pair<double, double> p;
  double px, dpx;
  int n = c.size();

  px = c[0] * std::pow(x, n - 1);
  for (int i = 1; i < n; ++i) {
    px = px + c[i] * std::pow(x, n - i - 1);
  }

  dpx = (n - 1) * c[0] * std::pow(x, n - 2);
  for (int i = 1; i < n - 1; ++i) {
    dpx = dpx + (n - i - 1) * c[i] * std::pow(x, n - i - 2);
  }
  
  p.first = px;
  p.second = dpx;
  
  return p;
}
/* SAM_LISTING_END_1 */
// END

/* SAM_LISTING_BEGIN_2 */
bool polyTestTime(unsigned int d) {
  
  bool ret = false;
  double x = .123;
  std::pair<double, double> p, p_naive;
  int repeats = 10;
  
  std::cout << std::setw(10) << "n" << std::setw(25) << "Horner scheme:" 
            << std::setw(25) << "Monomial approach:" << "\n"
            << " ================================================================\n";
  // TO DO (6-1.d): Compare the output of evaldp() and evaldp_naive() and
  // tabulate their runtimes for n=1,2,4,...,2^d.
  // START
  double TOL = 1E-8;
  int max_n = std::pow(2,d);
  VectorXd c = VectorXd::LinSpaced(max_n,1,max_n);
  for(int n = 1; n <= max_n; n*=2) {
    Timer tm_slow, tm_fast;
    
    for(int r=0; r<repeats; r++) {
      tm_slow.start();
      p_naive = evaldp_naive(c.head(n),x);
      tm_slow.stop();
        
      tm_fast.start();
      p = evaldp(c.head(n),x);
      tm_fast.stop();
    }
    
    double duration_slow = tm_slow.duration()/repeats;
    double duration_fast = tm_fast.duration()/repeats;
    
    std::cout << std::setw(10) << n << std::setw(25) << duration_fast 
              << std::setw(25) << duration_slow << "\n";
    
    // Check correctness
    double err_p = std::abs(p.first - p_naive.first);
    double err_dp = std::abs(p.second - p_naive.second);
    if ((err_p >= TOL) || (err_dp >= TOL)) {
      std::cout << "Discrepancy at n = " << n << std::endl;
      return false;
    }
  }
  ret = true;
  // END
  return ret;
}
/* SAM_LISTING_END_2 */

/*!
 * \brief dipoleval
 * \param t
 * \param y
 * \param x
 * \return
 */
/* SAM_LISTING_BEGIN_3 */
VectorXd dipoleval(const VectorXd &t, const VectorXd &y, const VectorXd &x) {
  assert(t.size() == y.size() && "t and y must have same size!");
  VectorXd ret; ret.resizeLike(x); // return value
  // TO DO (6-1.e): Use the Aitken-Neville algorithm to evaluate at x 
  // the derivative of the interpolating polynomial of the data (t,y).
  // START
  // Non-vectorized version: loop over all evaluation points
  for (int i = 0; i < x.size(); ++i) {
    VectorXd p(y);
    VectorXd dP = VectorXd::Zero(y.size());
    // Aitken-Neville outer loop
    for (int im = 1; im < y.size(); ++im) {
      for (int i0 = im - 1; i0 >= 0; --i0) {
        dP(i0) = (p(i0 + 1) + (x(i) - t(i0)) * dP(i0 + 1) - p(i0) -
                  (x(i) - t(im)) * dP(i0)) /
                 (t(im) - t(i0));
        p(i0) = ((x(i) - t(i0)) * p(i0 + 1) - (x(i) - t(im)) * p(i0)) /
                (t(im) - t(i0));
      }
    }
    ret(i) = dP(0);
  }
  // END
  return ret;
}
/* SAM_LISTING_END_3 */

/*!
 * \brief dipoleval_alt
 * \param t
 * \param y
 * \param x
 * \return
 */
/* SAM_LISTING_BEGIN_4 */
VectorXd dipoleval_alt(const VectorXd &t, const VectorXd &y,
                       const VectorXd &x) {
  assert(t.size() == y.size() && "t and y must have same size!");
  VectorXd ret;// return value
  // TO DO (6-1.f): Use polyfit() and polyval() to evaluate at x 
  // the derivative of the interpolating polynomial of the data (t,y).
  // START
  int n = y.size();
  VectorXd P = polyfit(t, y, n - 1).head(n - 1);

  for (int i = 0; i < P.size(); ++i) {
    P(i) *= n - 1 - i;
  }
  ret = polyval(P, x);
  // END
  return ret;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
bool testDipolEval(void) {
  bool ret = false;
  VectorXd dPx, dPx_alt;
  // TO DO (6-1.f): Compare dipoleval() and dipoleval_alt() for the data
  // given on the problem sheet.
  // START
  double TOL = 1E-8;
  int n = 10; // polynomial degree
  VectorXd x = VectorXd::LinSpaced(100,-1,4);
  VectorXd t = VectorXd::LinSpaced(n+1,0,3);
  VectorXd y = sin(t.array()).matrix();
  dPx = dipoleval(t,y,x);
  dPx_alt = dipoleval_alt(t,y,x);
  double err = (dPx-dPx_alt).norm();
  if (err < TOL) ret = true;
  // END
  return ret;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
void plotPolyInterpolant(const std::string &filename) {
  plt::figure();
  // TO DO (6-1.g): Plot the derivative of the polynomial that
  // interpolates the data from (6-1.f).
  // Hint: matplotlibcpp (here = plt) implements Python's
  // matplotlib functions such as figure(), plot(), xlabel(), title(),...
  // You can use main.cpp of the LinearDataFit problem as a reference.
  // START
  VectorXd dPx;
  int n = 10; // polynomial degree
  VectorXd x = VectorXd::LinSpaced(100,0,3);
  VectorXd t = VectorXd::LinSpaced(n+1,0,3);
  VectorXd y = sin(t.array()).matrix();
  dPx = dipoleval(t,y,x);
  plt::plot(x,dPx);
  plt::title("Derivative of interpolating polynomial");
  plt::xlabel("t");
  plt::ylabel("p'(t)");
  // END
  plt::savefig("./cx_out/" + filename + ".png");
}
/* SAM_LISTING_END_6 */
