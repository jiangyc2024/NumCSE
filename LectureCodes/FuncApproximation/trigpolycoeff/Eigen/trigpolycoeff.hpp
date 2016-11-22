# include <iostream>
# include <Eigen/Dense>
using Eigen::VectorXd;
using Eigen::MatrixXd;

/* SAM_LISTING_BEGIN_0 */
//Computes expansion coefficients of trigonometric polyonomials \eqref{eq:trigpreal}
// IN : \texttt{t} = vector of nodes \Blue{$t_{0},\ldots,t_n\in[0,1[$}
//      \texttt{y} =  vector of data \Blue{$y_{0},\ldots,y_{n}$}
// OUT: pair of vectors storing the basis expansion coefficients \Blue{$\alpha_j$}, \Blue{$\beta_j$}, see \cref{def:trip}
std::pair<VectorXd,VectorXd>
trigpolycoeff(const VectorXd& t, const VectorXd& y) {
  const unsigned N = y.size(), n = (N-1)/2;
  if (N % 2 == 0) throw "Number of points must be odd!\n";

  // build system matrix M
  MatrixXd M(N, N);
  M.col(0) = VectorXd::Ones(N);
  for (unsigned c = 1; c <= n; ++c) {
    M.col(c) = ( 2*M_PI*c*t ).array().cos().matrix();
    M.col(n + c) = ( 2*M_PI*c*t ).array().sin().matrix();
  }
  // solve LSE and extract coefficients \Blue{${\alpha_{j}}$} and \Blue{$\beta_{j}$}
  VectorXd c = M.lu().solve(y);
  return std::pair<VectorXd,VectorXd>(c.head(n + 1),c.tail(n));
}
/* SAM_LISTING_END_0 */


/* SAM_LISTING_BEGIN_1 */
//Computes expansion coefficients of trigonometric polyonomials \eqref{eq:trigpreal}
// IN : \texttt{t} = vector of nodes \Blue{$t_{0},\ldots,t_n\in[0,1[$}
//      \texttt{y} =  vector of data \Blue{$y_{0},\ldots,y_{n}$}
//      \texttt{a,b} are used to save the expansion coefficients \Blue{$\alpha_j$}, \Blue{$\beta_j$} 
void trigpolycoeff(const VectorXd& t, const VectorXd& y, VectorXd& a, VectorXd& b) {
  const unsigned N = y.size();
  if (N % 2 == 0) {
    std::cerr << "Number of points must be odd!\n";
    return;
  }
  const unsigned n = (N - 1) / 2;

  // build system matrix M
  MatrixXd M(N, N);
  M.col(0) = VectorXd::Ones(N);
  for (unsigned c = 1; c <= n; ++c) {
    M.col(c) = ( 2*M_PI*c*t ).array().cos().matrix();
    M.col(n + c) = ( 2*M_PI*c*t ).array().sin().matrix();
  }

  // solve LSE and extract coefficients \Blue{$\mathbf{\alpha}$} and \Blue{$\mathbf{\beta}$}
  VectorXd c = M.lu().solve(y);
  a = c.head(n + 1);
  b = c.tail(n);
}
/* SAM_LISTING_END_1 */
