# include <iostream>
# include <Eigen/Dense>
using Eigen::VectorXd;
using Eigen::MatrixXd;

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
  std::cout << M << "\n";

  // solve LSE and extract coefficients \Blue{$\mathbf{\alpha}$} and \Blue{$\mathbf{\beta}$}
  VectorXd c = M.lu().solve(y);
  a = c.head(n + 1);
  b = c.tail(n);
}
