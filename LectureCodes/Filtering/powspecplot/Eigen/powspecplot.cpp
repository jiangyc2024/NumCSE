# include <unsupported/Eigen/FFT>
# include <figure/figure.hpp>
# include "signalgen.hpp"
using Eigen::VectorXd;
using Eigen::VectorXcd;

int main() {
  VectorXd y = signalgen();
  Eigen::FFT<double> fft;
  VectorXcd c = fft.fwd(y);
  // imaginary part is 0, but we have to cut it off explicitly
  VectorXd p = c.cwiseProduct(c.conjugate()).real();

  mgl::Figure fig;
  fig.title("power spectrum");
  fig.bar(p.head(32));
  fig.xlabel("Index k of Fourier coefficient");
  fig.ylabel("|c_k|^2");
  fig.save("spec");
  return 0;
}
