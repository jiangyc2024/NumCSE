# include <Eigen/Dense>

void circul(Eigen::MatrixXd& C, const Eigen::VectorXd& a) {
  const long n = a.size();
  C = Eigen::MatrixXd::Zero(n, n);
  C.diagonal().array() += a(0);
  for (long i = 1; i < n; ++i) {
    C.diagonal(-i).array() += a(i);
    C.diagonal(n-i).array() += a(i);
  }
}
