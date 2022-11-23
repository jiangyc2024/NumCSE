# include <iostream>
# include "./blur.hpp"

int main() {
  const long n = 11, N = 7;
  MatrixXd P = MatrixXd::Ones(n, n),
           S = MatrixXd::Ones(N, N),
           C = blur(P, S);
  std::cout << "C:\n" << C << "\n";
  return 0;
}
