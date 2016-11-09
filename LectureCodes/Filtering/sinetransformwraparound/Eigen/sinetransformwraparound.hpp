#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

using Eigen::VectorXd;

/* SAM_LISTING_BEGIN_0 */
// Simple sine transform of \Blue{$\Vy\in\bbR^n-1$} into \Blue{$\Vc\in\bbR^{n-1}$} by \eqref{fft:sinft}
void sinetrfwrap(const VectorXd &y, VectorXd& c)
{
  VectorXd::Index n = y.size()+1;
  // Create wrapped vector \Blue{$\wt{\Vy}$}
  VectorXd yt(2*n); yt << 0,y,0,-y.reverse();
  
  Eigen::VectorXcd ct;
  Eigen::FFT<double> fft; // DFT helper class
  fft.SetFlag(Eigen::FFT<double>::Flag::Unscaled);
  fft.fwd(ct,yt);
  
  const std::complex<double> v(0,2); // factor \Blue{$2\iu$}
  c = (-ct.middleRows(1,n-1)/v).real();
}
/* SAM_LISTING_END_0 */
