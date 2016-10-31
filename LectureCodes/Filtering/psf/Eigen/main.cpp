# include <iostream>
# include "./psf.hpp"

int main() {
  MatrixXd S;
  const long L = 2;
  psf(L, S);
  std::cout << "S:\n" << S << "\n";
  return 0;
}
