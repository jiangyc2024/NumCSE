#include "cross.hpp"

using namespace Eigen;

int main(void) {
  double T = 1.;
  int N = 1;
  
  Vector3d y0;
  y0 << 0.1, 0.2, 0.4;
  
  auto f = [] (const Vector3d &y) -> Vector3d {
    Vector3d fy;
    fy << y(0)*y(1), y(1)*y(2), y(2)-y(0);
    return fy;
  };
  
  auto Jf = [] (const Vector3d &y) {
    Matrix3d J;
    J << y(1),y(0),0,0,y(2),y(1),-1,0,1;
    return J;
  };
  // test implicit midpoint
  std::vector<VectorXd> test_imp = solve_imp_mid(f, Jf, T, y0, N);
  std::cout << "Implicit midpoint:\n" << test_imp.back() 
            << std::endl<< std::endl;
  
  // test linear implicit midpoint
  std::vector<VectorXd> test_lin = solve_lin_mid(f, Jf, T, y0, N);
  std::cout << "Implicit linear midpoint:\n" << test_lin.back() 
            << std::endl<< std::endl;
  
  tab_crossprod();
}
