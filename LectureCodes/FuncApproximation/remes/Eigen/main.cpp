# include "./remes.hpp"

int main () {
  auto f = [](double x){ return std::sin(x); };
  auto df = [](double x){ return std::cos(x); };
  VectorXd c;
  remez(f, df, 0, 1, 5, 1e-5, c);

  return 0;
}
