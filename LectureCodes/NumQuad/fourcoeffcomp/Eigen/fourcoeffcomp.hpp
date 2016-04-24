# include <iostream>
# include <vector>
# include <complex>
# include <unsupported/Eigen/FFT>

template <class Function>
void fourcoeffcomp(std::vector<std::complex<double>>& y, Function& c, const unsigned m, const unsigned ovsmpl = 2) {
  // Compute the Fourier coefficients \Blue{$y_{-m}, \ldots, y_{m}$} of the function
  // \Blue{$c:[0,1[ \mapsto \bbC$} using an \Magenta{oversampling factor} \texttt{ovsmpl}.
  // \texttt{c} must be a handle to a function \texttt{\symbol{64}(t)}, e.g. a lambda function
  const unsigned N = (2*m + 1)*ovsmpl; // number of quadrature points
  const double h = 1./N; 

  // evaluate function in N points
  std::vector<std::complex<double>> c_eval(N);
  for (unsigned i = 0; i < N; ++i) {
    c_eval[i] = c(i*h);
  }
 
  // inverse discrete fourier transformation
  Eigen::FFT<double> fft;
  std::vector<std::complex<double>> z;
  fft.inv(z, c_eval);

  // Undo oversampling and wrapping of Fourier coefficient array
  // -> y contains same values as z but in a different order:
  // y = [z(N-m+1:N), z(1:m+1)]
  y = std::vector<std::complex<double>>();
  y.reserve(N);
  for (unsigned i = N - m; i < N; ++i) {
    y.push_back(z[i]);
  }
  for (unsigned j = 0; j < m + 1; ++j) {
    y.push_back(z[j]);
  }
}
