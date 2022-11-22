#include <Eigen/Dense>
#include <complex>

namespace foursumnaive {


using Eigen::VectorXd; using Eigen::VectorXcd;

//! Approximate direct evaluation of Fourier sum according to \eqref{eq:fourtrunc}
//  \texttt{signal} is a handle to a function providing the \Blue{$y_k$}
//  \texttt{M} specifies truncation of series according to \eqref{eq:fourtrunc}
//  \texttt{N} is the number of equidistant evaluation points for \Blue{$c$} in \Blue{$[0,1[$}
template <class Function>
VectorXcd foursumnaive(const Function& signal, int M, int N) {
  const std::complex<double> i(0,1); // imaginary unit
  // Evaluation points for Fourier sum \Blue{$c$}
  VectorXd t = VectorXd::LinSpaced(N,0,1-1./N); 
  VectorXcd c(N); c.setConstant(signal(0));
  const VectorXcd omega = (-2.*M_PI*i*t.array()).exp().matrix();
  VectorXcd omp = omega;
  VectorXcd omm = VectorXcd::Ones(N).cwiseQuotient(omega);

  // Inefficient direct summation of Fourier series
  for (int k = 1; k <= M; ++k) {
    c += signal(k)*omp + signal(-k)*omm;
    omp = omp.cwiseProduct(omega);          
    omm = omm.cwiseQuotient(omega);
  }
  return c;
}


} //namespace foursumnaive