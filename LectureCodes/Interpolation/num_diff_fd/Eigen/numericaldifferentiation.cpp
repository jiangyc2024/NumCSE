# include <vector>
# include <cmath>
# include <figure/figure.hpp>

# define SCALAR long double

template <typename Function>
void diff(const double& x, Function& f, Function& df) {
  std::vector<SCALAR> error, h;
  // build vector of difference quotient width
  for (int i = -1; i >= -61; i -= 5) {
    h.push_back( std::pow(2, i) );
  }
  for (int j = 0; j < h.size(); ++j) {
    const double df_approx = (f(x + h[j]) - f(x))/h[j];
    error.push_back( std::abs(df_approx - df(x)) );
  }

  mgl::Figure fig;
  fig.title("One sided-diff.quotient");
  fig.xlabel("h");
  fig.ylabel("error");
  fig.setlog(true, true);
  fig.plot(h, error, "+-");
  fig.save("expnumdiff.eps");
}

int main() {
  auto f = [](double x){ return std::exp(x); };
  auto df = f;
  diff(1.1, f, df);
  return 0;
}
