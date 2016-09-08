#include <Eigen/Dense>
#include <Eigen/QR>
using Eigen::VectorXd;
using Eigen::MatrixXd;

// Solution of linear least squares problem \eqref{eq:LSQ1} by means of QR-decomposition
// Note: \Blue{$\VA\in\bbR^{m,n}$} with \Blue{$m>n$}, \Blue{$\operatorname{rank}(\VA) = n$} is assumed
void qrlsqsolve(const MatrixXd& A, const VectorXd& b, VectorXd& x, double& res) {
  const unsigned m = A.rows(), n = A.cols();

  MatrixXd Ab(m, n + 1); Ab << A, b; // Ab = [A,b]
  
  // R = triu(qr(Ab,0)), QR-decomposition of extended matrix \label{qrl:1}
  MatrixXd R = Ab.householderQr().matrixQR().template triangularView<Eigen::Upper>(); 

  MatrixXd R_nn = R.block(0, 0, n, n); // R\_nn = R(1:n,1:n)
  // for the keyword "template", see http://eigen.tuxfamily.org/dox/TopicTemplateKeyword.html
  // R\_nn \ R(1:n,n+1), \Blue{$\wh{\Vx} = (\VR)_{1:n,1:n}^{-1}(\VQ^T\Vb)_{1:n}$}
  x = R_nn.template triangularView<Eigen::Upper>().solve(R.block(0, n, n, 1));     
  res = R(n, n); // \Blue{$= \N{\VA\wh{\Vx}-\Vb}_2$} (why ?) \label{qrl:2}
}
