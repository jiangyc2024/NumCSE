# include <iostream>
# include <Eigen/Dense>
# include <Eigen/Sparse>
# include <Eigen/SparseLU>

using Eigen::VectorXd;
using Eigen::SparseMatrix;

//**                                     **//
//**              Base class             **//
//**                                     **//

class Interpolation {
  public:
    Interpolation() {}
    virtual double operator()(const double x) const;
    virtual VectorXd operator()(const VectorXd& x) const;
    virtual ~Interpolation() {}
  protected:
    VectorXd t_, h_, y_, c_;
};

double Interpolation::operator()(const double x) const
{
  // find out between which nodes x is
  long i;
  for (i = 0; i < t_.size() - 1; ++i){
    if (t_(i) <= x && x <= t_(i + 1))
      break;
  }

  // evaluate
  const double tau = (x - t_(i))/h_(i);
  const double tau2 = tau*tau,
               tau3 = tau*tau*tau;
  const double res = y_(i)*(1 - 3*tau2 + 2*tau3) + y_(i + 1)*(3*tau2 - 2*tau3) + h_(i)*c_(i)*(tau - 2*tau2 + tau3) + h_(i)*c_(i + 1)*(tau3 - tau2);
  return res;
}

VectorXd Interpolation::operator()(const VectorXd& x) const {
  VectorXd y(x.size());
  for (long i = 0; i < x.size(); ++i) {
    y(i) = operator()(x(i));
  }
  return y;
}

//**                                     **//
//** Complete cubic spline interpolation **//
//**                                     **//

class CCSI : public Interpolation {
  public:
    CCSI (const VectorXd& t, const VectorXd& y, const double c0, const double cn);
    // operator() already implemented in base class Interpolation
    virtual ~CCSI() {};
};

CCSI::CCSI(const VectorXd& t, const VectorXd& y, const double c0, const double cn)
{
  y_ = y;
  t_ = t;
  const long n = t.size() - 1;
  h_ = (t.tail(n).array() - t.head(n).array()).matrix(); // size: n 
  VectorXd b = (1./h_.array()).matrix(); // size: n
  VectorXd a = 2*(b.head(n - 1).array() + b.tail(n - 1).array()).matrix(); // size: n - 1

  // build rhs
  Eigen::VectorXd rhs(n - 1);
  for (long i = 0; i < n - 1; ++i){
    rhs(i) = 3*( (y(i + 1) - y(i))/(h_(i)*h_(i)) + (y(i + 2) - y(i + 1))/(h_(i + 1)*h_(i + 1)) );
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

  Eigen::SparseLU<SparseMatrix<double>> solver;
  solver.compute(A);
  VectorXd cpart_ = solver.solve(rhs);
  c_ = Eigen::VectorXd(n + 1);
  c_ << c0, cpart_, cn;
} 


//**                                    **//
//** Natural cubic spline interpolation **//
//**                                    **//

class NCSI : public Interpolation {
  public:
    NCSI (const VectorXd& t, const VectorXd& y);
    // operator() already implemented in base class
    ~NCSI() {}
};

NCSI::NCSI(const VectorXd& t, const VectorXd& y)
{
  y_ = y; // accessing members of base class
  t_ = t;  
  const long n = t.size() - 1;
  h_ = (t.tail(n).array() - t.head(n).array()).matrix(); // size: n 
  VectorXd b = (1./h_.array()).matrix(); // size: n
  VectorXd a = 2*(b.head(n - 1).array() + b.tail(n - 1).array()).matrix(); // size: n - 1

  // build rhs -- size: n + 1
  Eigen::VectorXd rhs(n + 1);
  rhs(0) = 3*(y(1) - y(0))/(h_(1)*h_(1));
  for (long i = 0; i < n - 1; ++i){
    rhs(i + 1) = 3*( (y(i + 1) - y(i))/(h_(i)*h_(i)) + (y(i + 2) - y(i + 1))/(h_(i + 1)*h_(i + 1)) );
  }
  rhs(n) = 3*(y(n) - y(n - 1))/(h_(n - 1)*h_(n - 1));
  // modify according to complete cubic spline

  // build sparse system matrix -- size: (n + 1) x (n + 1)
  typedef Eigen::Triplet<double> T;
  std::vector<T> triplets;
  triplets.reserve(3*(n + 1));

  // first pair of entries initialized separately as the diagonal entry differs from all others
  triplets.push_back( T(0, 0, 2*b(0)) );
  triplets.push_back( T(0, 1, b(1)) );
  triplets.push_back( T(1, 0, b(0)) );

  // initalize middle part
  for (int i = 1; i < n; ++i) { // attention: index starts at 1
    triplets.push_back( T(i, i, a(i - 1)) ); // diagonal entries a_1 to a_(n - 1)
    triplets.push_back( T(i, i + 1, b(i)) ); // super-diagonal entries 
    triplets.push_back( T(i + 1, i, b(i - 1)) ); // sub-diagonal entries
  }

  // initialize last diagonal entry which is not covered by the main loop
  triplets.push_back( T(n, n, 2*b(n - 1)) );

  SparseMatrix<double> A(n + 1,n + 1);
  A.setFromTriplets(triplets.begin(), triplets.end());

  Eigen::SparseLU<SparseMatrix<double>> solver;
  solver.compute(A);
  c_ = solver.solve(rhs);
} 

//**                                     **//
//** Periodic cubic spline interpolation **//
//**                                     **//

class PCSI : public Interpolation {
  public:
    PCSI(const VectorXd& t, const VectorXd& y);
    ~PCSI() {}
};

PCSI::PCSI(const VectorXd& t, const VectorXd& y)
{
  y_ = y; // accessing members of base class
  t_ = t;  
  const long n = t.size() - 1;
  h_ = (t.tail(n).array() - t.head(n).array()).matrix(); // size: n 
  VectorXd b = (1./h_.array()).matrix(); // size: n
  VectorXd a = 2*(b.head(n - 1).array() + b.tail(n - 1).array()).matrix(); // size: n - 1

  // build rhs -- size: n 
  VectorXd rhs(n);
  for (long i = 0; i < n - 1; ++i){
    rhs(i) = 3*( (y(i + 1) - y(i))/(h_(i)*h_(i)) + (y(i + 2) - y(i + 1))/(h_(i + 1)*h_(i + 1)) );
  }
  rhs(n - 1) = 3*( (y(1) - y(0))/(h_(1)*h_(1)) + (y(n) - y(n - 1))/(h_(n - 1)*h_(n - 1)) );

  // build sparse system matrix -- size: n x n
  typedef Eigen::Triplet<double> T;
  std::vector<T> triplets;
  triplets.reserve(3*n);

  // initalize matrix - except last bottom right, bottom left and top bottom entries
  for (int i = 0; i < n - 1; ++i) { // attention: index starts at 1
    triplets.push_back( T(i, i, a(i)) ); // diagonal entries a_1 to a_(n - 1)
    triplets.push_back( T(i, i + 1, b(i + 1)) ); // super-diagonal entries 
    triplets.push_back( T(i + 1, i, b(i)) ); // sub-diagonal entries
  }

  // initialize last entries
  triplets.push_back( T(0, n - 1, b(0)) );
  triplets.push_back( T(n - 1, 0, b(0)) );
  triplets.push_back( T(n - 1, n - 1, 2*(b(0) + b(n - 1))) );

  SparseMatrix<double> A(n, n);
  A.setFromTriplets(triplets.begin(), triplets.end());

  std::cout << Eigen::MatrixXd(A) << "\n";

  Eigen::SparseLU<SparseMatrix<double>> solver;
  solver.compute(A);
  VectorXd c_part = solver.solve(rhs); // c_part = c1, ..., cn according to periodic splines: c0 = cn

  c_ = VectorXd(n + 1); // c_ = c0, ..., cn with c0 = cn
  c_ << c_part(n - 1), c_part;
}
