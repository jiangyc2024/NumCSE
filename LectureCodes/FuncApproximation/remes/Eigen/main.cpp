# include "./remes.hpp"

int main () {
  auto f = [](double x){ return 1/(1 + x*x); };
  auto df = [](double x){ return -2*x/( (1 + x*x)*(1 + x*x) ); };
  VectorXd c;
  remez(f, df, 0, 1, 5, 1e-10, c);
  std::cout << "RESULT ===============================================\n"
            << c.transpose() << "\n";

  return 0;
}
