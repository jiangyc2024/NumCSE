# include <cmath>
# include <vector>
# include <Eigen/Dense>
# include <figure/figure.hpp>
# include "gaussquad.hpp"

template <class Function>
void gaussConvCV(const Function& f, const double I_ex) {
  // transform integrand
  auto g = [f](double x) {
    return x*f(std::sin(x))*std::cos(x);
  };

  std::vector<double> evals, // used to save no. of quad nodes
                      error; // used to save the error
  const unsigned N = 50; // max. no. of nodes

  for (unsigned n = 1; n <= N; ++n) {
    QuadRule qr = gaussquad(n); // get quadrature rule
    // do not forget to transform quadratur rule!
    Eigen::VectorXd w = qr.weights * M_PI/2; 
    Eigen::VectorXd c = ( -M_PI/2 + M_PI/2*(qr.nodes.array() + 1) ).matrix();

    // evaluate g at quadrature nodes c
    Eigen::VectorXd gc = c.unaryExpr(g);
    // same as $I = \sum_{i=1}^n w_i g(c_i)$
    double I = w.dot(gc);  
    evals.push_back(n); // save no. of quadrature nodes
    error.push_back(std::abs(I - I_ex)); // save error

  }

  // create convergence plots
  mgl::Figure fig;
  fig.title("Gauss quadrature convergence");
  fig.setlog(false, true); // lin-log scaling
  fig.plot(evals, error, " +b").label("Error"); // plot error
  fig.xlabel("No. of quadrature nodes");
  fig.ylabel("|Error|");
  fig.legend();
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
