#include <iostream>
#include "dft2d.hpp"

using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
int main() {
  using Comp = complex<double>;
  const VectorXcd::Index m=7,n = 5;

  {
  // Test 2D convolution
    MatrixXd Y(m,n), X(m,n);
    MatrixXcd Z1(m,n),Z2(m,n);
    for (int k=0;k<m;k++)
      for (int j=0;j<n;j++) {
	Y(k,j) = std::min(k,j);
	X(k,j) = k+j;
      }
    pmconv_basic(X,Y,Z1);
    pmconv(X,Y,Z2);
    cout << "X = " << endl << X << endl;
    cout << "Y = " << endl << Y << endl;
    cout << "Z1 = " << endl << Z1 << endl;
    cout << "Z2 = " << endl << Z2 << endl;
    return 0;
  }
  
  {
  // Test: 2D DFT
  // Initialize matrix for testing 2D DFT
  // MATLAB: Y = triu((1:m)'*(1:n))+i*(ones(m,n));
  MatrixXcd Y = MatrixXcd::Constant(m,n,Comp(0.0,1.0));
  for (int k=0;k<m;k++)
    for (int j=0;j<n;j++)
      { if (k<=j) Y(k,j) += Comp(k*j,0.0); }

  MatrixXcd C(m,n);
  fft2(C,Y);
  cout << "Y = " << endl << Y << endl;
  cout << "C = " << endl << C << endl;
  ifft2(Y,C);
  cout << "Y = " << endl << Y << endl;
  return 0;
  }
}
/* SAM_LISTING_END_0 */
