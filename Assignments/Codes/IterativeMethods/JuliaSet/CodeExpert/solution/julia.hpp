#include <Eigen/Dense>
#include "matplotlibcpp.h"

using namespace Eigen;
namespace plt = matplotlibcpp;

/* SAM_LISTING_BEGIN_0 */
Vector2d juliaF(const Vector2d &z) {
  Vector2d Fz = Vector2d::Zero();
  // TODO: (9-8.c) Implement the function $F:R^2\to R^2$ such that F(z)=0 is
  // equivalent to z^3 = 1 (interpreting z as a complex number).
  // START
  Fz << z(0) * z(0) * z(0) - 3.0 * z(0) * z(1) * z(1) - 1.0,
      3.0 * z(0) * z(0) * z(1) - z(1) * z(1) * z(1);
  // END
  return Fz;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
Matrix2d juliaDF(const Vector2d &z) {
  Matrix2d DFz = Matrix2d::Zero();
  // TODO: (9-8.c) Implement the Jacobian of F at z.
  // START
  const double a = 3.0 * (z(0) * z(0) - z(1) * z(1));
  const double b = 6.0 * z(0) * z(1);
  DFz << a, -b, b, a;
  // END
  return DFz;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void julia(void) {
  // Exact solutions of $z^3 = 1$ for complex $z$.
  // Use 2-dimensional real vectors to represent complex numbers.
  Vector2d z1, z2, z3;
  z1 << 1, 0;
  z2 << -0.5, 0.5 * std::sqrt(3);
  z3 << -0.5, -0.5 * std::sqrt(3);

  // Tolerance and maximum number of iterations for Newton's method.
  constexpr double tol = 1e-4;
  constexpr unsigned int N_it = 20;

  // The image will have a resolution of res*res.
  constexpr unsigned int res = 780;
  // We imagine that we have a res*res grid Z on the square [-2,2] with
  // Z(0,0)=(-2,-2),    Z(0,res-1)=(-2,2),
  // Z(res-1,0)=(2,-2), Z(res-1,res-1)=(2,2), etc.
  // C(i,j) is the color assigned to the (i,j)-th point on the grid.
  MatrixXd C(res, res);

  // TODO: (9-8.d) Fill C with real numbers, such that each entry C(i,j) corresponds to
  // which of the roots (z1, z2, or z3) Newton's method converges
  // (if it converges), when using Z(i,j) as a starting point.
  // The values in C will be interpreted as colors below.
  // Additional: You may choose different shades of each color to indicate
  // how many iterations were needed to reach the tolerance.
  // Hint: To speed up runtimes, start with a low value for the parameter res.
  // START
  C.setZero();
  for (unsigned int i = 0; i < res; i++) {
    for (unsigned int j = 0; j < res; j++) {
      // Newton iteration starting at Z(i,j)=z=(x,y).
      const double x = -2 + 4.0 * i / (res - 1);
      const double y = -2 + 4.0 * j / (res - 1);
      Vector2d z(x, y);
      for (unsigned int k = 0; k < N_it; k++) {
        z -= juliaDF(z).lu().solve(juliaF(z));
        // Check if we are sufficiently close to any of the three roots:
        if ((z - z1).norm() < tol) {
          C(i, j) = 1.0 - 1.0 * k / N_it;
          break;
        }
        if ((z - z2).norm() < tol) {
          C(i, j) = 2.0 - 1.0 * k / N_it;
          break;
        }
        if ((z - z3).norm() < tol) {
          C(i, j) = 3.0 - 1.0 * k / N_it;
          break;
        }
      }
    }
  }
  // END
  // Axis labels
  std::vector<double> ticks(5);
  for (unsigned int i = 0; i < 5; ++i) {
    ticks[i] = i * (res - 1) / 4;
  }
  std::vector<std::string> labels{"-2", "-1", "0", "1", "2"};

  plt::figure();
  // Need to transpose C because plt::imshow() uses
  // first index for y axis and second index for x axis.
  // Different colormaps can be found at
  // https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
  plt::imshow(C.transpose(), {{"cmap", "viridis"}, {"origin", "lower"}});
  plt::colorbar();
  plt::xticks(ticks, labels);
  plt::yticks(ticks, labels);
  plt::title("Julia set for $z^3$ = 1");
  plt::savefig("./cx_out/julia.png");
}
/* SAM_LISTING_END_2 */
