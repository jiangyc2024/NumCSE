# include "deblur.hpp"
# include "../../blur/Eigen/blur.hpp"
# include "../../psf/Eigen/psf.hpp"
# include <unsupported/Eigen/KroneckerProduct>

int main() {
  MatrixXd M(3,3); M << 8,1,6,3,5,7,4,9,2;
  MatrixXd P = Eigen::kroneckerProduct(M, MatrixXd::Ones(2,2)),
           S; psf(1, S);
  MatrixXd C = blur(P, S);
  std::cout << "Original: \n" << P << "\n";
  std::cout << "Blurred: \n" << C << "\n";

  MatrixXd D = deblur(C,S).real();
  std::cout << "Deblurred: \n" << D << "\n";

  return 0;
}
