# include <cmath>
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "gaussquad.hpp"

template <class Function>
void gaussConvCV(const Function& f, const double I_ex) {
  // transform integrand
  auto g = [f](double x) {
    //
    // TODO: transform integrand as computed in previous exercise
    //
  };

  std::vector<double> evals, // used to save no. of quad nodes
                      error; // used to save the error
  const unsigned N = 50; // max. no. of nodes

  for (unsigned n = 1; n <= N; ++n) {
    //
    // TODO: get quadrature rule and compute the error
    //
  }

  // create convergence plots
  mgl::Figure fig;
  fig.title("Gauss quadrature convergence");
  //
  // TODO: create plot (use lin-log scaling!)
  //
  fig.save("GaussConvCV");
}

int main() {
  const double I_ex = 0.870267525725852642;
  // define $f(x) = \sinh x$
  std::function<double(double)> f = [](double x) {
    return std::sinh(x);
  };
  gaussConvCV(f, I_ex);
  return 0;
}
