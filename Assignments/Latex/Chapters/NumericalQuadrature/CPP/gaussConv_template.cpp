# include <cmath>
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "gaussquad.hpp"

// compute the integral I = \int arcsin(x)f(x)dx, x=-1,..,1
template <class Function>
void gaussConv(const Function& fh, const double I_ex) {

  // define integrand
  auto f = [fh](double x) {
    return std::asin(x)*fh(x);
  };

  std::vector<double> evals, // used to save no. of quad nodes
                      error; // used to save the error
  const unsigned N = 50; // max. no. of nodes
  
  for (unsigned n = 1; n <= N; ++n) {
    //
    // TODO: compute the error for n quadrature nodes
    //
  }

  // create convergence plots
  mgl::Figure fig;
  fig.title("Gauss quadrature convergence");
  //
  // TODO: create plots (use log-log scaling!)
  //
  fig.save("GaussConv");
}

int main() {
  // define f(x) = sinh(x)
  const double I_ex = 0.870267525725852642;
  std::function<double(double)> fh = [](double x) {
    return std::sinh(x);
  };
  gaussConv(fh, I_ex);
}
