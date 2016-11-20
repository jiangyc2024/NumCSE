# include <cmath>
# include <Eigen/Dense>
# include <unsupported/Eigen/FFT>
# include "stem.hpp" // simple classed which uses Figure for plotting
using Eigen::VectorXd; using Eigen::VectorXcd;

/* SAM_LISTING_BEGIN_0 */
// Visualize limit \Blue{$m\to\infty$} for a \Blue{$2m+!$}-periodic signal and
// its discrete Fourier transform ``squeezed'' into \Blue{$[0,1]$}.
int main() {
  // range of plot for visualization of discrete signal
  const int N = 3*257;
  // function defining discrete signal
  auto yfn = [](VectorXd k) {
    return (1/(1 + k.array()*k.array())).matrix();
  };
  // loop over different periods \Blue{$2^l+1$}
  for (int mpow = 4; mpow <= 7; ++mpow) {
    int m = std::pow(2, mpow);
    VectorXd ybas = yfn(VectorXd::LinSpaced(2*m+1, -m, m));
    const int Ncp = std::floor(N/(2*m+1));
    const int n = ybas.size();

    VectorXd y(n*Ncp);
    for (int i = 0; i < Ncp; ++i) 
      y.segment(n*i, n) = ybas;

    // Stem plots vertical lines from 0 to the y-value
    Stem s1;
    s1.title = "Period of signal y_i = " + std::to_string(2*m+1);
    s1.xlabel = "Index i of sampling instance";
    s1.ylabel = "y_i";
    s1.file = "persig" + std::to_string(mpow); // where to save
    s1.plot(VectorXd::LinSpaced(n*Ncp, -n*Ncp/2., n*Ncp/2.), y, "r");

    Eigen::FFT<double> fft;
    // DFT of wrapped signal (one period)
    VectorXd yconc(n); yconc << ybas.tail(m+1), ybas.head(m);
    VectorXcd c = fft.fwd(yconc);
    Stem s2;
    s2.title = "DFT of period " + std::to_string(2*m+1) + " signal";
    s2.xlabel = "t_k";
    s2.ylabel = "c(t_k)";
    s2.file = "persigdft" + std::to_string(mpow);
    s2.plot(VectorXd::LinSpaced(2*m+1, 0, 1), c.real(), "b");
  }
  return 0;
}
/* SAM_LISTING_END_0 */
