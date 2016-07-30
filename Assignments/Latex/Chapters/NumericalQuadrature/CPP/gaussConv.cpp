# include <cmath>
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "gaussQuad.hpp"

// computing the integral $\int_a^b f(x) dx \approx \sum_{i=1}^n w_i f(c_i)$
// for given quadrature rule $\{(w_i, x_i)\}_{i=1}^n$
template <class Function>
double integrate(const QuadRule& qr, const Function& f) {
  double I = 0;
  for (unsigned i = 0; i < qr.weights.size(); ++i) {
    I += qr.weights(i) * f(qr.nodes(i));
  }
  return I;
}

template <class Function>
void gaussConv(const Function& f, const double I_ex) {
  std::vector<double> evals, // used to save no. of quad nodes
                      error; // used to save the error
  const unsigned N = 50; // max. no. of nodes
  
  for (unsigned n = 1; n <= N; ++n) {
    QuadRule qr = gaussQuad(n); // get quadrature rule
    double I = integrate(qr, f); // compute integral
    evals.push_back(n); // save no. of quadrature nodes
    error.push_back(std::abs(I - I_ex)); // save error
  }

  // create convergence plots
  mgl::Figure fig;
  fig.title("Gauss quadrature convergence");
  fig.setlog(true, true); // log-log scaling
  fig.plot(evals, error, " +r").label("Error"); // plot error
  fig.fplot("x^(-3)", "k--").label("O(n^3)"); // reference line
  fig.xlabel("No. of quadrature nodes");
  fig.ylabel("|Error|");
  fig.legend();
  fig.save("GaussConv");
}

int main() {
  const double I_ex = 0.8702675180053843773; 
  // define $f(x) = \arcsin x \sinh x$
  std::function<double(double)> f = [](double x) {
    return std::asin(x)*std::sinh(x);
  };
  gaussConv(f, I_ex);
  return 0;
}
