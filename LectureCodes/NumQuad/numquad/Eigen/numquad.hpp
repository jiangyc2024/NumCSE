# include <vector>
# include <Eigen/Dense>
# include "polyfit.hpp"
# include "gaussquad.hpp"

template <class Function>
std::vector<double> numquad(Function& F, const double a, const double b, const unsigned N, const std::string mode) {
  // saving the results in this vector
  std::vector<double> res;
  res.reserve(N);

  if (mode == "gauss") {
    // using the struct QuadRule from gaussquad.hpp

    for (unsigned n = 1; n <= N; ++n) {
      QuadRule qr = gaussquad(n); // get nodes and weights
      // gaussquad returns nodes on [-1, 1] we must transform them to [a, b]
      Eigen::VectorXd x = (a + (b - a)/2*(qr.nodes.array() + 1)).matrix();

      // evaluate function at given nodes
      Eigen::VectorXd y(n);
      for (unsigned j = 0; j < n; ++j) {
        y(j) = F(x(j));
      }

      // compute integral using built-in dot function from Eigen
      // the (b - a)/2 term comes from transforming the weights from
      // [-1, 1] to [a, b] (sum of weights == length of interval)
      res.push_back( (b - a)/2*y.dot(qr.weights) );
    }
    return res;
  }

  // weights for clenshaw curtis quadrature
  Eigen::VectorXd w(N + 1);
  for (int i = N + 1; i >= 1; --i) {
    w(N + 1 - i) = (std::pow(b, i) - std::pow(a, i))/i;
  }
  
  for (unsigned n = 1; n <= N; ++n) {
    Eigen::VectorXd x(n + 1);
    // choose chebychev nodes on [a, b]
    if (mode == "chebychev") {
      for (unsigned j = 0; j <= n; ++j) {
        x(j) =  a + (b - a)/2*( std::cos( (2.*j + 1)/(2.*n + 2)*M_PI ) + 1 );
      }
    }
    else { // choose equidistant nodes on [a, b]
      x = Eigen::VectorXd::LinSpaced(n + 1, a, b);
    }

    // evaluate function at nodes
    Eigen::VectorXd y(n + 1);
    for (unsigned j = 0; j <= n; ++j) {
      y(j) = F(x(j));
    }

    // get monomial basis coefficients using polyfit
    Eigen::VectorXd coeffs = polyfit(x, y, n);

    // compute integral
    res.push_back( coeffs.dot(w.tail(n + 1)) );
  }
  return res;
}
