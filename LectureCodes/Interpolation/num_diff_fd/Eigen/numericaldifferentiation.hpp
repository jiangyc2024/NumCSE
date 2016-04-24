# include <vector>
# include <cmath>
# include <string>
# include <figure/figure.hpp>

# define SCALAR long double

// numerical differentiation using difference quotient \Blue{$f'(x) = \lim_{h\rightarrow 0} \frac{f(x + h) - f(x)}{h}
// IN: f = function to derive, 
//     df = exact derivative (to compare results),
//     name = string of function name (for plot filename)
// OUT: plot of error will be saved as "<name>numdiff.eps"

template <class Function, class Derivative>
void diff(const double x, Function& f, Derivative& df, const std::string name) {
  std::vector<SCALAR> error, h;
  // build vector of difference quotient width
  for (int i = -1; i >= -61; i -= 5) {
    h.push_back( std::pow(2, i) );
  }
  for (int j = 0; j < h.size(); ++j) {
    // compute approximate solution using difference quotient
    const double df_approx = (f(x + h[j]) - f(x))/h[j],
    // compute relative error
                 rel_err = std::abs( (df_approx - df(x))/df(x) );
    std::cout << "h = " << h[j] << "  " << rel_err << "\n";
    error.push_back(rel_err);
  }

  mgl::Figure fig;
  fig.title("One sided-diff.quotient");
  fig.xlabel("h");  fig.ylabel("error");
  fig.setlog(true, true);
  fig.plot(h, error, "+-");
  fig.save(name + "numdiff.eps");
}
