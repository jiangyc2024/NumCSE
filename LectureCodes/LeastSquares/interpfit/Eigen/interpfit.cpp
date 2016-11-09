///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Xiaolin Guo, Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

# include <Eigen/Dense>
# include <figure/figure.hpp>
# include <polyfit.hpp> // NCSE's polyfit (equivalent to Matlab's)
# include <polyval.hpp> // NCSE's polyval (equivalent to Matlab's)

// Comparison of polynomial interpolation and polynomial fitting
// (``Quick and dirty'', see \ref{subsec:interp-algorithms})
int main() {
  // use C++ lambda functions to define runge function \Blue{$f(x) = \frac{1}{1 + x^2}$}
  auto f = [](const Eigen::VectorXd& x){ 
    return ( 1./(1 + x.array()*x.array()) ).matrix();
  };

  const unsigned d = 10; // Polynomial degree
  Eigen::VectorXd tip(d + 1); // \Blue{$d + 1$} nodes for interpolation
  for (unsigned i = 0; i <= d; ++i)
    tip(i) = -5 + i*10./d;

  Eigen::VectorXd tft(3*d + 1); // \Blue{$3d + 1$} nodes for polynomial fitting
  for (unsigned i = 0; i <= 3*d; ++i) 
    tft(i) = -5 + i*10./(3*d);

  Eigen::VectorXd ftip = f(Eigen::VectorXd::Ones(2));
  Eigen::VectorXd pip = polyfit(tip, f(tip), d),  // Interpolating polynomial (deg = \Blue{$d$})
                  pft = polyfit(tft, f(tft), d); // Fitting polynomial (deg = \Blue{$d$})

  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(1000,-5,5);
  mgl::Figure fig;
  fig.plot(x, f(x), "g|").label("Function f");
  fig.plot(x, polyval(pip, x), "b").label("Interpolating polynomial");
  fig.plot(x, polyval(pft, x), "r").label("Fitting polynomial");
  fig.plot(tip, f(tip), " b*");
  fig.save("interpfit");


  return 0;
}
