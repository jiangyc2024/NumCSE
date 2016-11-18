# include <iostream>
# include <complex>
# include <Eigen/Dense>

// Local includes
# include "intpolyval_complex.hpp"
# include "ipvclass.hpp"

using Eigen::VectorXd;
using Eigen::VectorXcd;


/* SAM_LISTING_BEGIN_0 */
// Evaluation of trigonometric interpolant at numerous points
// IN : \texttt{t} = vector of nodes \Blue{$t_{0},\ldots,t_n\in[0,1[$}
//      \texttt{y} = vector of data \Blue{$y_{0},\ldots,y_{n}$}
//      \texttt{x} = vector of evaluation points \Blue{$x_{1},\ldots,x_{N}$}
// OUT : vector of values of the interpolant at points in \texttt{x}
VectorXd trigpolyval(const VectorXd& t, const VectorXd& y, const VectorXd& x) {
  using idx_t = VectorXd::Index;   
  using comp = std::complex<double>; 
  const idx_t N = y.size(); // Number of data points
  if (N % 2 == 0) throw std::runtime_error("Number of points must be odd!");
  const idx_t n = (N - 1)/2;
  const std::complex<double> M_I(0,1); // imaginary unit
  // interpolation nodes and evalutation points on unit circle
  VectorXcd tc = ( 2*M_PI*M_I*t ).array().exp().matrix(),
            xc = ( 2*M_PI*M_I*x ).array().exp().matrix();
  // Rescaled values, according to \Blue{$q(t) = e^{-2\pi int}\cdot p(e^{2\pi it})$}, see \eqref{eq:trigpcomp}
  VectorXcd z = ((2*n*M_PI*M_I*t).array().exp() * y.array()).matrix();
  // Evaluation of interpolating polynomial on unit circle using the
  // barycentric interpolation formula in \Blue{$\bbC$}, see \cref{cpp:ipvclass}
  BarycPolyInterp<comp> Interpol(tc);
  VectorXcd p = Interpol.eval<VectorXcd>(z,xc);
  // Undo the scaling, see \eqref{eq:scale} 
  VectorXcd qc = ((-2*n*M_PI*M_I*x).array().exp() * p.array()).matrix();
  return (qc.real()); // imaginary part is zero, cut it off
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
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
/* SAM_LISTING_END_1 */
