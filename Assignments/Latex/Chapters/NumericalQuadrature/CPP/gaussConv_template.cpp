# include <cmath>
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "gaussQuad.hpp"

template <class Function>
void gaussConv(const Function& f, const double I_ex) {
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
  fig.save("conv");
}

int main() {
  // define f(x) = asin(x)*sinh(x)
  const double I_ex = 0.8702675180053843773; 
  std::function<double(double)> f = [](double x) {
    return std::asin(x)*std::sinh(x);
  };
  gaussConv(f, I_ex);
}
