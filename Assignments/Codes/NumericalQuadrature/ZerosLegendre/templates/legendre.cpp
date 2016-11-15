#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

// Evaluate the Legendre polynomials and its derivatives in vector $x$ 
// using the 3-term recursion formulae.
// The outputs are the matrices $Lx$ and $DLx$.
/* SAM_LISTING_BEGIN_0 */
void legvals(const VectorXd& x, MatrixXd& Lx, MatrixXd& DLx) {
  const int n = Lx.cols()-1;
  const int N = x.size();
    // TODO: evaluate Legendre polynomials
}
/* SAM_LISTING_END_0 */

// Evaluate $P_n(x)$ for a scalar $x$ and integer $n$.
/* SAM_LISTING_BEGIN_1 */
double Pnx(const double x, const int n) {
  VectorXd Px(n+1);
    // TODO: evaluate $P_n(x)$
  return Px(n);
}
/* SAM_LISTING_END_1 */

// Find the Gauss points using the secant method without regula falsi.
/* SAM_LISTING_BEGIN_2 */
MatrixXd gaussPts(const int n, const double rtol=1e-10,
							   const double atol=1e-12) {
  MatrixXd zeros(n,n);
    // TODO: secant method without regula falsi
  return zeros;
}
/* SAM_LISTING_END_2 */

// Find the Gauss points using the secant method with regula falsi.
// The standard secant method may be obtained
// by commenting out lines 92 and 93.
/* SAM_LISTING_BEGIN_3 */
MatrixXd gaussPts_regulaFalsi(const int n, const double rtol=1e-10,
										   const double atol=1e-12) {
  MatrixXd zeros(n,n);
    // TODO: secant method with regula falsi
  return zeros;
}
/* SAM_LISTING_END_3 */

int main() {
  const int n = 8;
  
  // Secant method without regula falsi
  std::cout << "---> Secant method without regula falsi\n";
  
  MatrixXd zeros = gaussPts(n);
  std::cout << "Zeros:\n" << zeros << "\n";
    
  for (int k = 1; k < n+1; ++k) {
    VectorXd xi = zeros.block(0, k-1, k, 1);
    MatrixXd Lx(k,n+1), DLx(k,n+1);
    legvals(xi, Lx, DLx);
    std::cout << "Values of the " << k
			  << "-th polynomial in the calculated zeros:\n"
              << Lx.col(k).transpose() << "\n";
  }
  
  // Secant method with regula falsi
  std::cout << "---> Secant method with regula falsi\n";
  
  MatrixXd zeros = gaussPts_regulaFalsi(n);
  std::cout << "Zeros:\n" << zeros << "\n";
    
  for (int k = 1; k < n+1; ++k) {
    VectorXd xi = zeros.block(0, k-1, k, 1);
    MatrixXd Lx(k,n+1), DLx(k,n+1);
    legvals(xi, Lx, DLx);
    std::cout << "Values of the " << k
			  << "-th polynomial in the calculated zeros:\n"
              << Lx.col(k).transpose() << "\n";
  }
}
