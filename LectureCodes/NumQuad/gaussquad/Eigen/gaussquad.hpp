# include <iostream>
# include <Eigen/Dense>
# include <Eigen/Eigenvalues>


struct QuadRule {
  Eigen::VectorXd nodes, weights;
};

void gaussquad (const unsigned n, QuadRule& qr)
{
  // Eigen::MatrixXd M(n,n); would yield a Matrix with nonzero entries!
  // therefore we must explicitly set them to zero
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
  for (unsigned i = 1; i < n; ++i){
    const double b = i/std::sqrt(4.*i*i - 1.);
    M(i, i - 1) = b;
    // line 19 is optional as the EV-Solver only references 
    // the lower triangular part of M
    // M(i - 1, i) = b; 
  }
  
  // using EigenSolver for symmetric matrices, exploiting the structure
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(M);

  qr.nodes = eig.eigenvalues();
  qr.weights = 2*eig.eigenvectors().topRows<1>().array().pow(2);
}
