# include <iostream>
# include <cmath>
# include <Eigen/Dense>
using Eigen::VectorXd;
using Eigen::MatrixXd;

void vander(const VectorXd& x, MatrixXd& V);

template <class Function>
VectorXd feval(const Function& f, const VectorXd& x);

// IN : f = handle to the Function, point evaluation
//      df = handle to the Derivative of f, point evaluation
//      a, b = interval boundaries 
//      d = degree of polynomial
// OUT: c = coefficients to interpolant in monomial basis

template <class Function, class Derivative>
void remez(const Function& f, const Derivative& df, const double a, const double b, const unsigned d, const double tol, VectorXd& c) {
  const unsigned n = 8*d; // number of sampling points
  VectorXd xtab = VectorXd::LinSpaced(n - 1, a, b); // points of sampling grid
  VectorXd ftab = feval(f, xtab); // function values at sampling grid
  double fsupn = ftab.cwiseAbs().maxCoeff(); // approximate supremum norm of \Blue{$f$}
  VectorXd dftab = feval(df, xtab); // derivative values at sampling grid

  // The vector xe stores the current guess for the alternants
  // initial guess is Chebychev alternant \eqref{remez:chebalt}
  const double h = M_PI/(d + 1);
  VectorXd xe(d + 2);
  for (int i = 0; i < d + 2; ++i) {
    xe(i) = (a + b)/2 + (a - b)/2*std::cos(h*i);
  }

  VectorXd fxe = feval(f, xe); // f evaluated at alternants

  const unsigned maxit = 10;
  // Main iteration loop of remez algorithm
  for (unsigned k = 0; k < maxit; ++k) {
    // Interpolation at \Blue{$d+2$} points xe with deviations \Blue{$\pm\delta$}
    // Algorithm uses monomial basis, which is \textit{not} optimal
    MatrixXd V, A;
    vander(xe, V); // get Vandermonde matrix for nodes xe
    for (unsigned i = 0; i < V.rows(); ++i) {
      V(i, d + 1) = std::pow(-1, i); // modify last column according to algorithm
    }

    // solve for coefficients and error
    VectorXd c = V.lu().solve(fxe),
             cd(d);
    for (unsigned i = 0; i < d; ++i) {
      cd(i) = (d - i)*c(d - i);
    }
    cd(
    // find extrema
    // (i) Find initial guesses for the inner extremes by sampling
    //     track sign changes of the derivative of the approximation error
    

  }

}

template <class Function, class Derivative>
std::vector<double> findExtrema(const Function& f, const Derivative& df);

// [1, x_1, ..., x_1^n-1]
//    .      .         .
//    .        .       .
// [1, x_n, ..., x_n^n-1]
void vander(const VectorXd& x, MatrixXd& V) {
  const unsigned n = x.size();
  V = MatrixXd::Ones(n, n);
  for (unsigned c = 1; c < n; ++c) {
    V.col(c) = x.array().pow(c).matrix();
  }
}

template <class Function>
VectorXd feval(const Function& f, const VectorXd& x) {
  VectorXd fx(x.size());
  for (unsigned i = 0; i < x.size(); ++i) {
    fx(i) = f(x(i));
  }
  return fx;
}
