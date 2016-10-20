#include <iostream>
#include "dft2d.hpp"

using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
int main() {
  using Comp = complex<double>;
  const VectorXcd::Index m=7,n = 5;
  // Initialize matrix for testing
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
  
}
/* SAM_LISTING_END_0 */
