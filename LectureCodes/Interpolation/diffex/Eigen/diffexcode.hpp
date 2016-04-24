# include <vector>
# include <Eigen/Dense>

using Eigen::VectorXd;

template <class Function>
// \texttt{f}: handle of a function defined in a neighbourhood of \Blue{$x \in \mathbb{R}$}
// \texttt{x}: point at which approximate derivative is desired
// \texttt{h0}: initial distance from \texttt{x}
// \texttt{rtol}: relative target tolerance, \texttt{atol}: absolute tolerance
double diffex(Function& f, const double x, const double h0, const double rtol, const double atol) {
  const unsigned nit = 10;
  VectorXd h(nit), y(nit);
  h(0) = h0;
  y(0) = (f(x + h0) - f(x - h0))/(2*h0);

  // using Aitken-Neville scheme with \Blue{$x=0$}, see Code~\ref{AitkenNeville} 
  for (unsigned i = 1; i < nit; ++i) {
    // create interpolation data
    h(i) = h(i - 1)/2;
    y(i) = ( f(x + h(i)) - f(x - h(i)) )/h(i - 1);

    // aitken neville update
    for (int k = i - 1; k >= 0; --k) { // careful: dont use type(k) = unsigned!
      y(k) = y(k + 1) - (y(k + 1) - y(k))*h(i) / (h(i) - h(k));
    }
    
    // termination of extrapolation when desired tolerance is achieved
    const double errest = std::abs(y(1) - y(0)); // error indicator
    if ( errest < rtol*std::abs(y(0)) || errest < atol ){ // \label{de:1}
      break;
    }
  }
  return y(1);
}
