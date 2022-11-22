#include <Eigen/Dense>
#include <cmath>
#include <unsupported/Eigen/FFT>

namespace foursum {


using Eigen::VectorXcd;
using Eigen::VectorXd;
using std::ceil;

/* SAM_LISTING_BEGIN_0 */
// evaluate scalar function with a vector
//  DFT based approximate evaluation of Fourier series
//  \texttt{signal} is a functor providing the \Blue{$y_k$}
//  \texttt{M} specifies truncation of series according to \eqref{eq:fourtrunc}
//  \texttt{N} is the number of equidistant evaluation points for \Blue{$c$} in
//  \Blue{$[0,1[$}.
template <class Function>
VectorXcd foursum(const Function &signal, int M, int N) {
  typedef Eigen::Index index_t;
  const int m = 2 * M + 1; // length of the signal
  // sample signal
  // VectorXd y = feval(signal, VectorXd::LinSpaced(m, -M, M));
  VectorXd y = VectorXd::LinSpaced(m, -M, M).unaryExpr(signal);
  // Ensure that there are more sampling points than terms in series
  const int l = static_cast<int>(m > N ? ceil(static_cast<double>(m) / N) : 1);
  N *= l;
  // \Hyperlink{ZEROPAD}{\cor{Zero padding}} and wrapping of signal, see
  // \cref{cpp:freqfilter}
  VectorXd y_ext = VectorXd::Zero(N);
  y_ext.head(M + 1) = y.tail(M + 1);
  y_ext.tail(M) = y.head(M);
  // Perform DFT and decimate output vector
  Eigen::FFT<double> fft;
  const Eigen::VectorXcd k = fft.fwd(y_ext);
  Eigen::VectorXcd c(N / l);
  for (int i = 0; i < N / l; ++i) {
    c(i) = k(static_cast<index_t>(i) * static_cast<index_t>(l));
  }
  return c;
}
/* SAM_LISTING_END_0 */


} //namespace foursum