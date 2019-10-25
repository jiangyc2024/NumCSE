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

/* SAM_LISTING_BEGIN_0 */
// START

// END
/* SAM_LISTING_END_0 */



/*
 * @brief Evaluate a polynomial and its derivative using a naive implementation
 * @param[in] vector c of size $n$, coefficients of the polynomial p
 * @param[in] double x, where the polynomial has to be evaluated
 * @param[out] pair containing p(x),p'(x)
 */
// TO DO (6-1.b): Implement the function evaldp_naive().
/* SAM_LISTING_BEGIN_1 */
// START

// END
/* SAM_LISTING_END_1 */


/* SAM_LISTING_BEGIN_2 */
bool polyTestTime(unsigned int d) {
  
  bool ret = false;
  double x = .123;
  std::pair<double, double> p, p_naive;
  int repeats = 10;
  
  std::cout << std::setw(10) << "n" << std::setw(25) << "Horner scheme:" 
            << std::setw(25) << "Monomial approach:" << "\n"
            << " ================================================================\n";
  // TO DO (6-1.d): Compare the output of evaldp() and evaldp_naive() for 
  // x = .123 and tabulate their runtimes for n=2,4,...,2^d.
  // START
  
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
  
  // END
  plt::savefig("./cx_out/" + filename + ".png");
}
/* SAM_LISTING_END_6 */
