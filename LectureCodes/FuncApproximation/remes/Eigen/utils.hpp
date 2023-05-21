// \emph{utils.hpp}

# include "polyval.hpp" // from NumCSE/Utils
# include <Eigen/Dense>
# include <algorithm>
# include <vector>

namespace remez {


using Eigen::VectorXd;
using Eigen::MatrixXd;

// compute the Vander matrix to the vector x:
// [$x_1^{n - 1}$ ...      $x_1^0$]
// [ .           .          . ]
// [$x_n^{n - 1}$ ...      $x_n^0$]
inline MatrixXd vander(const VectorXd& x) {
  const Eigen::Index n = x.size();
  MatrixXd V(n,n);
  for (int c = 0; c < n; ++c) {
    V.col(n - 1 - c) = x.array().pow(c).matrix();
  }
  return V;
}

// returns f(x) = [f(x1), .., f(xn)] as Eigen::VectorXd
template <class Function>
VectorXd feval(const Function& f, const VectorXd& x) {
  VectorXd fx(x.size());
  for (unsigned i = 0; i < x.size(); ++i) {
    fx(i) = f(x(i));
  }
  return fx;
}

// returns all negative entries in x as Eigen::VectorXd
inline VectorXd findNegative(const VectorXd& x) {
  std::vector<double> negs;
  for (int i = 0; i < x.size(); ++i) {
    if(x(i) < 0) {
      negs.push_back(i);
    }
  }
  VectorXd res = Eigen::Map<VectorXd>(negs.data(), static_cast<Eigen::Index>(negs.size()));
  return res;
}

// select from x all elements with index in ind, equivalent to 
// Matlab's x(ind), where ind is a vector
inline VectorXd select(const VectorXd& x, const VectorXd& ind) {
  VectorXd res(ind.size());
  for (int i = 0; i < ind.size(); ++i) {
    res(i) = x( static_cast<Eigen::Index>(ind(i)));
  }
  return res;
}

// returns how the indices have changed after x has been sorted,
// equivalent to Matlab's [xSorted, ind] = sort(x)
inline VectorXd sort_indices(const VectorXd& x) {
  VectorXd ind = VectorXd::LinSpaced(x.size(), 0, static_cast<double>(x.size() - 1));
  std::sort(ind.begin(), ind.end(),
      [&x](int i1, int i2) { return x(i1) < x(i2); });
  return ind;
}


} //namespace remez
