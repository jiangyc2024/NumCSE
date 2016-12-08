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
  for (int j = 0; j < N; ++j) {
    Lx(j,0) = 1.;
    Lx(j,1) = x(j);
    DLx(j,0) = 0;
    DLx(j,1) = 1.;
    for (int k = 2; k < n+1; ++k) {
      Lx(j,k) = (2*k-1.)/k*x(j)*Lx(j,k-1) - (k-1.)/k*Lx(j,k-2);
      DLx(j,k) = (2*k-1.)/k*Lx(j,k-1) + (2*k-1.)/k*x(j)*DLx(j,k-1) -
				 (k-1.)/k*DLx(j,k-2);
    }
  }
}
/* SAM_LISTING_END_0 */

// Evaluate $P_n(x)$ for a scalar $x$ and integer $n$.
/* SAM_LISTING_BEGIN_1 */
double Pnx(const double x, const int n) {
  VectorXd Px(n+1);
  Px(0) = 1.; Px(1) = x;
  for (int k=2; k<n+1; k++){
    Px(k) = (2*k-1.)/k*x*Px(k-1)-(k-1.)/k*Px(k-2);
  }
  return Px(n);
}
/* SAM_LISTING_END_1 */

// Find the Gauss points using the secant method without regula falsi.
/* SAM_LISTING_BEGIN_2 */
MatrixXd gaussPts(const int n, const double rtol=1e-10,
							   const double atol=1e-12) {
  MatrixXd zeros(n,n);
  double x0, x1, f0, f1, s;
  for (int k = 1; k < n+1; ++k) {
    for (int j = 1; j < k+1; ++j) {
      // Initial guesses
      if (j == 1) x0 = -1.;
      else      x0 = zeros(j-2, k-2);
      if (j == k) x1 = 1.;
      else      x1 = zeros(j-1, k-2);
            
      // Secant method
      f0 = Pnx(x0, k);
      for (int i = 0; i < 1e4; ++i) {
        f1 = Pnx(x1, k);
        s = f1*(x1-x0)/(f1-f0);
        x0 = x1; f0 = f1;
        x1 -= s;
        if ( (std::abs(s) < std::max(atol, rtol*std::min(std::abs(x0),
              std::abs(x1)))) ) {
          zeros(j-1, k-1) = x1;
          break;
        }
      }
    }
  }
  return zeros;
}
/* SAM_LISTING_END_2 */

// Find the Gauss points using the secant method with regula falsi.
// The standard secant method may be obtained
// by commenting out line 106.
/* SAM_LISTING_BEGIN_3 */
MatrixXd gaussPts_regulaFalsi(const int n, const double rtol=1e-10,
										   const double atol=1e-12) {
  MatrixXd zeros(n,n);
  double x0, x1, f0, f1, s;
  for (int k = 1; k < n+1; ++k) {
    for (int j = 1; j < k+1; ++j) {
      // Initial guesses
      if (j == 1) x0 = -1;
      else      x0 = zeros(j-2, k-2);
      if (j == k) x1 = +1;
      else      x1 = zeros(j-1, k-2);

      // Secant method
      f0 = Pnx(x0, k);
      for (int i = 0; i < 1e4; ++i) {
        f1 = Pnx(x1, k);
        s = f1*(x1-x0)/(f1-f0);
        if (Pnx(x1-s, k)*f1 < 0) // NEW LINE: regula falsi
          x0 = x1; f0 = f1;
        x1 -= s;
        if ( (std::abs(s) < std::max(atol, rtol*std::min(std::abs(x0),
			  std::abs(x1)))) ) {
          zeros(j-1, k-1) = x1;
          break;
        }
      }
    }
  }
  return zeros;
}
/* SAM_LISTING_END_3 */

int main() {
  const int n = 8;
  MatrixXd zeros;
  
  // Secant method without regula falsi
  std::cout << "---> Secant method without regula falsi\n";
  
  zeros = gaussPts(n);
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
  
  zeros = gaussPts_regulaFalsi(n);
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
