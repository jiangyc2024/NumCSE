# ifndef GRID_HPP
# define GRID_HPP

# include <Eigen/Dense>

typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
std::pair<Mat, Mat> meshgrid(Vec& a, Vec& b)
{
  long m = a.size();
  long n = b.size();
  Mat X(n,m), Y(n,m);
  for (int i = 0; i < n; ++i)
    X.block(i, 0, 1, m) = a.transpose();
  for (int j = 0; j < m; ++j)
    Y.block(0, j, n, 1) = b;
  return std::pair<Mat,Mat>(X,Y);
}

# endif
