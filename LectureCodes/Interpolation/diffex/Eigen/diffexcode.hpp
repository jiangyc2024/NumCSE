///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Chapter on data interpolation
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <vector>
# include <Eigen/Dense>

using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
// Extrapolation based numerical differentation
// with a posteriori error control
// \texttt{f}: handle of a function defined in a neighbourhood of \Blue{$x \in \mathbb{R}$}
// \texttt{x}: point at which approximate derivative is desired
// \texttt{h0}: initial distance from \texttt{x}
// \texttt{rtol}: relative target tolerance, \texttt{atol}: absolute tolerance
template <class Function>
double diffex(Function& f, const double x, const double h0,
	      const double rtol, const double atol) {
  const unsigned nit = 10; // Maximum depth of extrapolation
  VectorXd h(nit); h[0] = h0; // Widths of difference quotients
  VectorXd y(nit); // Approximations returned by difference quotients
  y[0] = (f(x + h0) - f(x - h0))/(2*h0); // Widest difference quotients

  // using \com{Aitken-Neville scheme} with \Blue{$x=0$}, see Code~\ref{AitkenNeville} 
  for (unsigned i = 1; i < nit; ++i) {
    // create data points for extrapolation
    h[i] = h[i-1]/2; // Next width half a big
    y[i] = ( f(x + h[i]) - f(x - h[i]) )/h(i - 1);
    // Aitken-Neville update
    for (int k = i - 1; k >= 0; --k) 
      y[k] = y[k+1] - (y[k+1]-y[k])*h[i]/(h[i]-h[k]);
    // termination of extrapolation when desired tolerance is reached
    const double errest = std::abs(y[1]-y[0]); // \com{error indicator}
    if ( errest < rtol*std::abs(y[0]) || errest < atol ) // \label{de:1}
      break;
  }
  return y[0]; // Return value extrapolated from largest number of difference quotients 
}
/* SAM_LISTING_END_0 */
