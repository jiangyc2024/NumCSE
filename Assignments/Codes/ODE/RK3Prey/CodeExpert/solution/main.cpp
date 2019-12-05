
#include "rkintegrator.hpp"
#include "rk3prey.hpp"

using namespace Eigen;

int main(void) {
  MatrixXd A(2,2);
  A << 0, 0, 1, 0;
  VectorXd b(2);
  b << 0.5, 0.5;
  auto f = [] (VectorXd y) {
    VectorXd fy(2);
    fy << -0.5*y(0), y(0)*y(1);
    return fy;
  };
  VectorXd y0(2);
  y0 << -1,1;
  
  RKIntegrator<VectorXd> Solver(A,b);
  auto Y = Solver.solve(f,2,y0,10);
  std::cout << "Test of RKIntegrator.solve() gives y(2) = " << Y[10].transpose() << std::endl;
  
  std::cout << "\nResults of RK3prey():\n";
  std::cout << "Convergence rate: " << std::round(std::abs(RK3prey())) << std::endl;
}
