# include "pconvfft.hpp"

VectorXcd myconv(const VectorXcd& h, const VectorXcd& x) {
  const long n = h.size();
  // Zero padding, cf. \eqref{eq:zeropad}
  VectorXcd hp(2*n - 1), xp(2*n - 1);
  hp << h, VectorXcd::Zero(n - 1);
  xp << x, VectorXcd::Zero(n - 1);
  // Periodic discrete convolution of length \Blue{$2n-1$}
  return pconvfft(hp, xp);
}
