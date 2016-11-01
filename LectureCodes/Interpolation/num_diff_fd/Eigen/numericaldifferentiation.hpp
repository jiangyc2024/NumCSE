///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): J. Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

# include <vector>
# include <cmath>
# include <string>
# include <figure/figure.hpp>

/* SAM_LISTING_BEGIN_0 */
// numerical differentiation using difference quotient
// \Blue{$f'(x) = \lim_{h\rightarrow 0} \frac{f(x + h) - f(x)}{h}$}
// IN: f (function object) = function to derive 
//     df = exact derivative (to compute error),
//     name = string of function name (for plot filename)
// OUT: plot of error will be saved as "<name>numdiff.eps"
template <class Function, class Derivative>
void diff(const double x, Function& f, Derivative& df, const std::string name) {
  std::vector<long double> error, h; 
  // build vector of widths of difference quotients
  for (int i = -1; i >= -61; i -= 5) h.push_back(std::pow(2,i));
  for (unsigned j = 0; j < h.size(); ++j) {
    // compute approximate solution using difference quotient
    double df_approx = (f(x + h[j]) - f(x))/h[j];
    // compute relative error
    double rel_err = std::abs((df_approx - df(x))/df(x));
    error.push_back(rel_err);
  }
/* SAM_LISTING_END_0 */
  mgl::Figure fig;
  fig.title("One sided-diff.quotient");
  fig.xlabel("h");  fig.ylabel("error");
  fig.setlog(true, true);
  fig.plot(h, error, "+-");
  fig.save(name + "numdiff.eps");
}
