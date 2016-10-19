# include <Eigen/Dense>
using Eigen::VectorXcd;

VectorXcd pconv(const VectorXcd& u, const VectorXcd& x) {
  const long n = x.size(); 
  VectorXcd z = VectorXcd::Zero(n);
  for (long k = 0; k < n; ++k) {
    for (long j = 0; j < n; ++j) {
      long ind = (k - j < 0 ? n + k - j : k - j);
      z(k) += u(ind)*x(j);
    }
  }
  return z;
}
