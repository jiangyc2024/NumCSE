# include <cmath>
# include <Eigen/Dense>
# include <Eigen/Eigenvalues>

struct QuadRule {
  Eigen::VectorXd nodes, weights;
};

QuadRule gaussquad(const unsigned n) {
  QuadRule qr;
  qr.nodes.resize(n);
  qr.weights.resize(n);
  if (n == 0) {
    qr.nodes(0) = 0;
    qr.weights(0) = 2;
  } else {
    Eigen::VectorXd b(n-1);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n,n);

    for(unsigned i = 1; i < n; ++i) {
      double d = i / std::sqrt(4.*i*i - 1.);
      J(i,i-1) = d;
      J(i-1,i) = d;
    }
    Eigen::EigenSolver<Eigen::MatrixXd> eig(J);

    qr.nodes = eig.eigenvalues().real();
    qr.weights = 2 * eig.eigenvectors().real().topRows<1>().cwiseProduct(eig.eigenvectors().real().topRows<1>());
  }
  return qr;
}
