# include <iostream>
# include <complex>
# include <Eigen/Dense>
# include "intpolyval_complex.hpp"

using Eigen::VectorXd;
using Eigen::VectorXcd;

/* SAM_LISTING_BEGIN_0 */
// Evaluation of trigonometric interpolant at numerous points
// IN : \texttt{t} = vector of nodes \Blue{$t_{0},\ldots,t_n\in[0,1[$}
//      \texttt{y} = vector of data \Blue{$y_{0},\ldots,y_{n}$}
//      \texttt{x} = vector of evaluation points \Blue{$x_{1},\ldots,x_{N}$}
//      \texttt{q} is used to save the value of the interpolant at \texttt{x}
void trigpolyval(const VectorXd& t, const VectorXd& y, const VectorXd& x, VectorXd& q) {
  const int N = y.size();
  if (N % 2 == 0) {
    std::cerr << "Number of points must be odd!\n";
    return;
  }
  const int n = (N - 1)/2;
  const std::complex<double> i(0,1); // imaginary unit
  // interpolation nodes and evalutation points on unit circle
  VectorXcd tc = ( 2*M_PI*i*t ).array().exp().matrix(),
            xc = ( 2*M_PI*i*x ).array().exp().matrix();
  // Rescaled values, according to \Blue{$q(t) = e^{-2\pi int}\cdot p(e^{2\pi it})$}
  VectorXcd z = ((2*n*M_PI*i*t).array().exp() * y.array()).matrix();
  // Evaluation of interpolating polynomial on unit circle, see Code~\ref{barycentricformula}
  VectorXcd p; intpolyval(tc, z, xc, p);
  // Undo the scaling, see \eqref{eq:scale} 
  VectorXcd qc = ((-2*n*M_PI*i*x).array().exp() * p.array()).matrix();
  q = qc.real(); // imaginary part is zero, cut it off
}
/* SAM_LISTING_END_0 */
