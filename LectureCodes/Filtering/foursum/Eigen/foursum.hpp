# include "feval.hpp" // evaluate scalar function with a vector
# include <cmath>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
using Eigen::VectorXd; using Eigen::VectorXcd;

//! Approximate evaluation of Fourier series
//  \texttt{signal} is a handle to a function providing the \Blue{$y_k$}
//  \texttt{M} specifies truncation of series according to \eqref{eq:fourtrunc}
//  // \texttt{N} is the number of equidistant evaluation points for \Blue{$c$} in \Blue{$[0,1[$}. 
template <class Function>
VectorXcd foursum(const Function& signal, int M, int N) {
  const int m = 2*M+1; // length of the signal
  // sample signal
  VectorXd y = feval(signal, VectorXd::LinSpaced(m, -M, M));
  // Ensure that there are more sampling points than terms in series
  int l;
  if (m > N) {
    l = std::ceil(double(m)/N); N *= l;
  }
  else
    l = 1;
  
  // \Hyperlink{ZEROPAD}{Zero padding} and wrapping of signal, see Code~\ref{mc:freqfilter}
  VectorXd y_ext = VectorXd::Zero(N);
  y_ext.head(M+1) = y.tail(M+1);
  y_ext.tail(M) = y.head(M);

  // Perform DFT and decimate output vector
  Eigen::FFT<double> fft;
  Eigen::VectorXcd k = fft.fwd(y_ext), c(N/l);
  for (int i = 0; i < N/l; ++i) {
    c(i) = k(i*l);
  }
  return c;
}
