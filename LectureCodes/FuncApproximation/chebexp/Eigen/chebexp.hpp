///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
# include <cassert>
# include <cmath>
# include <vector>
# include <array>
# include <complex> // needed to use complex std::vectors!
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT> // contains Eigen::FFT

using Eigen::VectorXd;
using Eigen::VectorXcd;

/* SAM_LISTING_BEGIN_0 */
// efficiently compute coefficients \Blue{$\alpha_j$} in the Chebychev expansion
// \Blue{$p = \sum\limits_{j=0}^{n} \alpha_j T_j$} of \Blue{$p\in\Cp_n$} based on values \Blue{$y_k$},
// \Blue{$k=0,\ldots,n$}, in Chebychev nodes \Blue{$t_k$}, \Blue{$k=0,\ldots,n$}
// IN:  values \Blue{$y_k$} passed in \texttt{y}
// OUT: coefficients \Blue{$\alpha_j$}
VectorXd chebexp(const VectorXd& y) {
  const int n = y.size() - 1;          // degree of polynomial
  const std::complex<double> M_I(0, 1); // imaginary unit
  // create vector \Blue{$\Vz$}, see \eqref{eq:tpipcext}
  VectorXcd b(2*(n + 1));
  const std::complex<double> om = -M_I*(M_PI*n)/((double)(n+1)); 
  for (int j = 0; j <= n; ++j) {
    b(j) = std::exp(om*double(j))*y(j); // this cast to double is necessary!!
    b(2*n+1-j) = std::exp(om*double(2*n+1-j))*y(j);
  }
  
  // Solve linear system \eqref{eq:tpiplse} with effort \com{$O(n\log n)$}
  Eigen::FFT<double> fft;  // \eigen's helper class for DFT
  VectorXcd c = fft.inv(b); // -> c = ifft(z), inverse fourier transform 
  // recover \Blue{$\beta_j$}, see \eqref{eq:tpiplse}
  VectorXd beta(c.size());
  const std::complex<double> sc = M_PI_2/(n + 1)*M_I;
  for (unsigned j = 0; j < c.size(); ++j) 
    beta(j) = ( std::exp(sc*double(-n+j))*c[j] ).real();
  // recover \Blue{$\alpha_j$}, see \eqref{eq:tpdft}
  VectorXd alpha = 2*beta.segment(n,n); alpha(0) = beta(n);
  return alpha;
}
/* SAM_LISTING_END_0 */

// Conversion from Chebychev basis to monomial basis
// Argument a provides coefficients in Chebychev basis
// returns vector of monomial coefficients 
Eigen::VectorXd chebexpToMonom(const Eigen::VectorXd &a) {
  const int n = a.size() -1; // degree of polynomial
  assert(n >= 0); 
  // vector for returning monomial coefficients of polynomials
  Eigen::VectorXd r{ Eigen::VectorXd::Zero(n+1) };
  r[0] = a[0];
  if (n > 0) {
    r[1] = a[1];
    if (n > 1) {
      // monomial coefficients of Chebychev polynomials 
      std::array<Eigen::VectorXd, 2> c { Eigen::VectorXd::Zero(n+1), Eigen::VectorXd::Zero(n+1)};
      c[0][0] = 1.0; // Lowest-degree Chebychev polynomial is constant = 1
      c[1][1] = 1.0; // First-degree Chebychev polynomial is t -> t
      for (int j=2; j<=n; ++j) {
	c[j%2][0] *= -1.0;
	for (int k=1; k <= j; ++k) 
	  c[j%2][k] = 2*c[(j+1)%2][k-1] - c[j%2][k];
        r += a[j]*c[j%2];
      }
    }
  }
  return r;
}
