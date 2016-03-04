# include <cmath>
# include "./numquad.hpp"

int main () {
  auto f = [](double x){ return x*x; };
  const double a = -1,
               b = 1;
  numquad(f, a, b, 10, "equidistant");
  numquad(f, a, b, 10, "chebychev");
  numquad(f, a, b, 10, "gauss");

  return 0;
}
