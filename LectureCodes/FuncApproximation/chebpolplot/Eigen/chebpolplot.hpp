///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "chebpolmult.hpp"

/* SAM_LISTING_BEGIN_0 */
// plots Chebychev polynomials up to degree \texttt{nmax} on \Blue{$[-1,1]$}
void chebpolplot(const unsigned nmax) {
  Eigen::RowVectorXd x = Eigen::RowVectorXd::LinSpaced(500, -1, 1); // evaluation points
  Eigen::MatrixXd V; chebpolmult(nmax, x, V); // get values of cheb. polynomials 

  mgl::Figure fig;
  // iterate over rows of \texttt{V}, which contain the values of the cheb. polynomials
  for (unsigned r = 0; r < nmax + 1; ++r) 
    fig.plot(x, V.row(r)).label("n = " + std::to_string(r));

  fig.title("First " + std::to_string(nmax) + " cheb. polynomials");
  fig.xlabel("x"); fig.ylabel("T_n(x)");
  fig.ranges(-1.1, 1.1, -1.1, 1.1);
  fig.legend(); fig.save("chebpols.eps");
}
/* SAM_LISTING_END_0 */

