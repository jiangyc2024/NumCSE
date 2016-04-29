# include <cmath>
# include "./numquad.hpp"

int main () {
  //auto f = [](double x){ return x; };
  //auto f = [](double x){ return 1/(1 + 25*x*x); };
  auto f = [](double x){ return std::sqrt(x); };
  const double a = 0,
               b = 1;
  //numquad(f, a, b, 20, "equidistant");
  //numquad(f, a, b, 20, "chebychev");
  numquad(f, a, b, 20, "gauss");

  return 0;
}
