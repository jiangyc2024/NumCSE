# include <iostream>
# include "trigipequid.hpp" // FuncApproximation/trigipequid/Eigen
# include <complex>
# include <Eigen/Dense>
using Eigen::VectorXcd; 
using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
// Evaluation of trigonometric interpolation polynomial through \Blue{$(\frac{j}{2n+1},y_j)$}, \Blue{$j=0,\ldots,2n$}
// in equidistant points \Blue{$\frac{k}{N}$}, \Blue{$k=0,N-1$}
// IN : \texttt{y} = vector of values to be interpolated
//      \texttt{q} (COMPLEX!) will be used to save the return values
void trigpolyvalequid(const VectorXd y, const int M, VectorXd& q) {
  const int N = y.size();
  if (N % 2 == 0) {
    std::cerr << "Number of points must be odd!\n";
    return;
  }
  const int n = (N - 1)/2;
  // computing coefficient \Blue{$\gamma_j$}, see \eqref{tip:FM}
  VectorXcd a, b;
  trigipequid(y, a, b);

  std::complex<double> i(0,1);
  VectorXcd gamma(2*n + 1);
  gamma(n) = a(0);
  for (int k = 0; k < n; ++k) {
    gamma(k) = 0.5*( a(n - k) + i*b(n - k - 1) );
    gamma(n + k + 1) = 0.5*( a(k + 1) - i*b(k) );
  }  

  // zero padding
  VectorXcd ch(M); ch << gamma, VectorXcd::Zero(M - (2*n + 1));

  // build conjugate fourier matrix
  Eigen::FFT<double> fft;
  const VectorXcd chCon = ch.conjugate(); 
  const VectorXcd v = fft.fwd(chCon).conjugate();

  // multiplicate with conjugate fourier matrix
  VectorXcd q_complex = VectorXcd(M);
  for (int k = 0; k < M; ++k) {
    q_complex(k) = v(k) * std::exp( -2.*k*n*M_PI/M*i );
  }
  // complex part is zero up to machine precision, cut off!
  q = q_complex.real();
}
/* SAM_LISTING_END_0 */
