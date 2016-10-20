#include <complex>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

using namespace std;
using namespace Eigen;

/* SAM_LISTING_BEGIN_0 */
int main() {
  using Comp = complex<double>;
  const VectorXcd::Index n = 5;
  VectorXcd y(n),c(n),x(n);
  y << Comp(1,0),Comp(2,1),Comp(3,2),Comp(4,3),Comp(5,4);
  FFT<double> fft; // DFT transform object
  c = fft.fwd(y);  // DTF of y, see \cref{def:DFT} 
  x = fft.inv(c);  // inverse DFT of c, see \eqref{eq:invdft}

  cout << "y = " << y.transpose() << endl
       << "c = " << c.transpose() << endl
       << "x = " << x.transpose() << endl;
  return 0;
}
/* SAM_LISTING_BEGIN_0 */
