# include "shapepresintp.hpp"

int main() {
  VectorXd t0 = VectorXd::LinSpaced(13, 0, 12),
           y0(13);
  y0 << 0, 0.3, 0.5, 0.2, 0.6, 1.2, 1.3, 1, 1, 1, 1, 0, -1;

  shapepresintp(t0, y0);
  //shapepresintp(u, z);
  return 0;
}
