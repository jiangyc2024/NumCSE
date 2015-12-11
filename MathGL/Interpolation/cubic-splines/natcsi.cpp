# include "natcsi.hpp"
# include <Eigen/Dense>
# include <Eigen/Sparse> 
# include <Eigen/SparseCholesky>
# include <vector>
# include <algorithm>
# include <cassert>

// PRE: sorted vector of t and corresponding interpolation values y
// POST: create object of NatCSI
NatCSI::NatCSI(const Eigen::VectorXd& t, const Eigen::VectorXd& y)
: t_(t), y_(y) 
{
  assert(t.size() == y.size());
  const std::size_t n = t.size() - 1; // t = t0, ..., tn -> #t = n+1
  // build helper-definitions
  const Eigen::VectorXd h = t.tail(n) - t.head(n); // #h = n 
  const Eigen::VectorXd b = (1./h.array()).matrix(); // #b = n
  const Eigen::VectorXd a = 2*(b.head(n-1) + b.tail(n-1));
  const Eigen::VectorXd diff_y = y.tail(n) - y.head(n);
  Eigen::VectorXd rhs_constr = (diff_y.array()/h.array()/h.array()).matrix();
  
  // build rhs
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);
  rhs.head(n) = rhs_constr;
  // need to save it temporarily otherwise Eigen starts overwriting the old vector while we still need it
  Eigen::VectorXd rhs_tail_tmp = rhs.tail(n); 
  rhs.tail(n) = rhs_tail_tmp + rhs_constr;

  // build system matrix
  typedef Eigen::Triplet<double> T;
  std::vector<T> triplets;
  triplets.reserve(3*(n + 1) - 2);

  // first row:
  triplets.push_back( T(0, 0, 2*b(0)) );
  triplets.push_back( T(0, 1, b(0)) );

  // rows 2 to n:
  for (std::size_t i = 1; i < n; ++i){
    triplets.push_back( T(i, i, a(i - 1)) );
    triplets.push_back( T(i, i - 1, b(i - 1)) );
    triplets.push_back( T(i, i + 1, b(i)) );
  }

  // last row:
  triplets.push_back( T(n, n - 1, b(n - 1)) );
  triplets.push_back( T(n, n, b(n - 1)) );

  Eigen::SparseMatrix<double> sysmat(n + 1, n + 1);
  sysmat.setFromTriplets(triplets.begin(), triplets.end());

  // solve LSE
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.analyzePattern(sysmat);
  solver.factorize(sysmat);
  c_ = solver.solve(rhs);
}

// PRE: x in between t_0 and t_end
double NatCSI::operator() (double x) const{
  // find out between which nodes x is
  unsigned int index;
  // the case of x being equal to the last node must be handled separately because the find_if
  // will fail due to the "<"
  if (x == t_(t_.size() - 1)){
    return y_(t_.size() - 1);
  }
  else
  {
    auto f = [x](double node){ return x < node; };
    auto index_pointer = std::find_if(t_.data(), t_.data() + t_.size(), f);
    index = index_pointer - t_.data();
  }
  double h = t_(index) - t_(index - 1);
  x = (x - t_(index - 1))/h;
  double a1 = y_(index) - y_(index - 1);
  double a2 = a1 - h*c_(index - 1);
  double a3 = h*c_(index) - a1 - a2;
  double s =  y_(index - 1) + (a1 + (a2 + a3*x)*(x - 1))*x;
  return s;
}
