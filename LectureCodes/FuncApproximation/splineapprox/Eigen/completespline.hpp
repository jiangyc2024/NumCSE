# include <Eigen/Dense>
# include <Eigen/Sparse>
# include <Eigen/SparseLU>

using Eigen::VectorXd;
using Eigen::SparseMatrix;


// perform complete cubic spline interpolation of data \Blue{$(t, y)$} and evaluate at points \texttt{x}
// implementation at the bottom
VectorXd spline(const VectorXd& t, const VectorXd& y, const double c0, const double cn, const VectorXd& x);

// get coefficients for complete cubic spline interpolation
// IN : \Blue{$(t, y)$} = set of interpolation data
//      \texttt{c0, cn} = slope at start and end point
// OUT: coefficients for complete cubic spline interpolation
VectorXd coeffsCSI(const VectorXd& t, const VectorXd& y, const double c0, const double cn)
{
  const long n = t.size() - 1;
  const VectorXd h = (t.tail(n).array() - t.head(n).array()).matrix(); // size: n 
  const VectorXd b = (1./h.array()).matrix(); // size: n
  const VectorXd a = 2*(b.head(n - 1).array() + b.tail(n - 1).array()).matrix(); // size: n - 1

  // build rhs
  Eigen::VectorXd rhs(n - 1);
  for (long i = 0; i < n - 1; ++i){
    rhs(i) = 3*( (y(i + 1) - y(i))/(h(i)*h(i)) + (y(i + 2) - y(i + 1))/(h(i + 1)*h(i + 1)) );
  }
  // modify according to complete cubic spline
  rhs(0) -= b(0)*c0;
  rhs(n - 2) -= b(n - 1)*cn;

  // build sparse system matrix
  typedef Eigen::Triplet<double> T;
  std::vector<T> triplets;
  triplets.reserve(3*n);
  for (int i = 0; i < n - 2; ++i){
    triplets.push_back( T(i, i, a(i)) );
    triplets.push_back( T(i, i + 1, b(i + 1)) );
    triplets.push_back( T(i + 1, i, b(i)) );
  }
  triplets.push_back( T(n - 2, n - 2, a(n - 2)) );
  SparseMatrix<double> A(n - 1,n - 1);
  A.setFromTriplets(triplets.begin(), triplets.end());

  // solve LSE and apply conditions
  Eigen::SparseLU<SparseMatrix<double>> solver;
  solver.compute(A);
  const VectorXd cpart = solver.solve(rhs);
  VectorXd c = Eigen::VectorXd(n + 1);
  c << c0, cpart, cn;

  return c;
}

// evaluate cubic spline at a point x
// IN : \texttt{x} = evaluation point
//      \Blue{$(t, y)$} = set of interpolation data
//      \texttt{c} = vector of coefficients
// OUT: interpolated value at \texttŧ{x}
double eval(const double x, const VectorXd& t, const VectorXd& y, const VectorXd& c) {
  // find out between which nodes x is
  long i, n = t.size() - 1;
  for (i = 0; i < n; ++i){
    if (t(i) <= x && x <= t(i + 1))
      break;
  }
  if (i == n) { // didnt find x
    --i; // but as we access h(i), i must be smaller than n!
  }

  // vector of meshwidths
  const VectorXd h = t.tail(n) - t.head(n);

  // evaluate
  const double tau = (x - t(i))/h(i);
  const double tau2 = tau*tau,
               tau3 = tau*tau*tau;
  const double res = y(i)*(1 - 3*tau2 + 2*tau3) + y(i + 1)*(3*tau2 - 2*tau3) + h(i)*c(i)*(tau - 2*tau2 + tau3) + h(i)*c(i + 1)*(tau3 - tau2);
  return res;
}

// evaluate cubic spline at a vector x
// IN : \texttt{x} = evaluation points
//      \Blue{$(t, y)$} = set of interpolation data
//      \texttt{c} = vector of coefficients
// OUT: interpolated values at point x
VectorXd eval(const VectorXd& x, const VectorXd& t, const VectorXd& y, const VectorXd& c) {
  VectorXd fx(x.size());
  for (long i = 0; i < x.size(); ++i) {
    fx(i) = eval(x(i), t, y, c);
  }
  return fx;
}

VectorXd spline(const VectorXd& t, const VectorXd& y, const double c0, const double cn, const VectorXd& x) {
  const VectorXd c = coeffsCSI(t, y, c0, cn);
  return eval(x, t, y, c);
  return c;
}
