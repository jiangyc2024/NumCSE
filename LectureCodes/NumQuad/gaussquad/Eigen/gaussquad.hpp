# include <Eigen/Dense>
# include <Eigen/Eigenvalues>


struct QuadRule {
  Eigen::VectorXd nodes;
  Eigen::VectorXd weights;
};

void gaussquad (const unsigned& n, QuadRule& qr)
{
  // Eigen::MatrixXd M(n,n); would yield a Matrix with nonzero entries!
  // therefore we must explicitly set them to zero
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
  for (unsigned i = 1; i < n; ++i){
    const double b = i/std::sqrt(4.*i*i - 1.);
    // set matrix entries of GW-Matrix
    M(i, i - 1) = b;
    M(i - 1, i) = b;
  }

  // intitialize EigenSolver with our Matrix
  Eigen::EigenSolver<Eigen::MatrixXd> eig(M);

  // the matrix is symmetric and has only real eigenvalues but as this is
  // a numeric procedure small complex parts can appear 
  // -> we need to cut them off
  qr.nodes = eig.eigenvalues().real();

  // same goes for the eigenvectors
  const Eigen::VectorXd temp = eig.eigenvectors().topRows<1>().real();
  qr.weights = 2*temp.cwiseProduct(temp);

  return;
}
