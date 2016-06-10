# include <cmath>
# include <vector>
# include <complex>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>

using Eigen::VectorXd;

// efficiently compute coefficients \Blue{$\alpha_j$} in the Chebychev expansion
// \Blue{$p = \sum\limits_{j=0}^{n} \alpha_j T_j$} of \Blue{$p\in\Cp_n$} based on values \Blue{$y_k$},
// \Blue{$k=0,\ldots,n$}, in Chebychev nodes \Blue{$t_k$}, \Blue{$k=0,\ldots,n$}
// IN:  values \Blue{$y_k$} passed in \texttt{y}
// OUT: coefficients \Blue{$\alpha_j$}
VectorXd chebexp(const VectorXd& y) {
  const int n = y.size() - 1; // degree of polynomial
  const std::complex<double> i(0, 1); // imaginary unit
  // create vector \Blue{$\Vz$} 
  std::vector<std::complex<double>> z(2*(n + 1));
  std::complex<double> k = -i*(M_PI*n)/(n+1.); // just a constant, no meaning
  for (int j = 0; j <= n; ++j) {
    z[j] = std::exp(k*double(j))*y(j); // this cast to double is necessary!!
    z[2*n+1-j] = std::exp(k*double(2*n+1-j))*y(j);
  }
  
  // Solve linear system \eqref{eq:tpiplse} with effort \Blue{$O(n\log n)$}
  std::vector<std::complex<double>> c; // vector to save result of ifft(z)
  Eigen::FFT<double> fft; 
  fft.inv(c, z); // -> c = ifft(z), inverse fourier transform 
  // recover \Blue{$\beta_j$}, see \eqref{eq:tpiplse}
  VectorXd b(c.size());
  k = M_PI_2/(n + 1)*i;
  for (int j = 0; j < c.size(); ++j) {
    b(j) = ( std::exp(k*double(-n+j))*c[j] ).real();
  }

  // recover \Blue{$\alpha_j$}, see \eqref{eq:tpdft}
  VectorXd a = 2*b.tail(n);
  a(0) = b(n);
  return a;
}
