# include <vector>
# include <cmath>
# include <string>
# include <figure/figure.hpp>

# define SCALAR long double

template <class Function, class Derivative>
void diff(const double& x, Function& f, Derivative& df, const std::string& name) {
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
