// **********************************************************************
// Eigen codes for testing dense linear algebra routines
// **********************************************************************

#include <Eigen/Dense>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <span>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

using std::cout;
using std::endl;

void lsesolve(Eigen::Index n = 7)
{
  // Initialize a special invertible matrices
  MatrixXd mat = MatrixXd::Identity(n,n) +
    VectorXd::Constant(n,1.0)*RowVectorXd::Constant(n,1.0);
  cout << "Matrix mat = " << endl << mat << endl;
  const MatrixXd utm = mat.triangularView<Eigen::Upper>();
  cout << "Matrix utm = " << endl << utm << endl;
  // Multiple right hand side vector
  const MatrixXd B = MatrixXd::Random(n,2);
  // Solve linear system using various decompositions
  const MatrixXd X = mat.lu().solve(B);
  const MatrixXd X2 = mat.fullPivLu().solve(B);
  const MatrixXd X3 = mat.householderQr().solve(B);
  const MatrixXd X4 = mat.llt().solve(B);
  const MatrixXd X5 = mat.ldlt().solve(B);
  cout << "|X2-X| = " << (X2-X).norm() << endl;
  cout << "|X3-X| = " << (X3-X).norm() << endl;
  cout << "|X4-X| = " << (X4-X).norm() << endl;
  cout << "|X5-X| = " << (X5-X).norm() << endl;
}

void reuselu(double tol = 1E-3, Eigen::Index n = 7)
{
  // Initialize a special invertible matrices
  const MatrixXd A = MatrixXd::Identity(n,n) +
    VectorXd::Constant(n,1.0)*RowVectorXd::Constant(n,1.0);
  // Request LU-decomposition
  auto A_lu_dec = A.lu();
  // Initial guess for inverse power iteration
  VectorXd xo = VectorXd::Zero(n);
  VectorXd xn = VectorXd::Random(n); 
  xn /= xn.norm();
  while ((xo-xn).norm() > xn.norm()*tol) {
    xo = xn;
    xn = A_lu_dec.solve(xo);
    xn /= xn.norm();
    cout << "x = [" << xn.transpose() << "]" << endl;
  }
}

template<class VecType, class MatType>
VecType invpowit(const Eigen::MatrixBase<MatType> &A,double tol)
{
  using index_t = typename MatType::Index;
  const index_t n = A.cols();
  const index_t m = A.rows();
  eigen_assert(n == m);
  // Request LU-decomposition
  auto A_lu_dec = A.lu();
  // Initial guess for inverse power iteration
  VecType xo = VecType::Zero(n);
  VecType xn = VecType::Random(n); 
  xn /= xn.norm();
  while ((xo-xn).norm() > xn.norm()*tol) {
    xo = xn;
    xn = A_lu_dec.solve(xo);
    xn /= xn.norm();
  }
  return(xn);
}

void invpowitdriver(double tol = 1E-3, Eigen::Index n = 7)
{
  // Initialize a special invertible matrices
  const MatrixXd A = MatrixXd::Identity(n,n) +
    VectorXd::Constant(n,1.0)*RowVectorXd::Constant(n,1.0);
  auto ev = invpowit<VectorXd>(A,tol);
  cout << "ev = [" << ev.transpose() << "]" << endl;
  VectorXd av = A*ev;
  cout << "ev/av = [" << (ev.array()/av.array()).transpose() << "]" << endl;
}

void fn( Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> &A,
	const Eigen::VectorXd &v,Eigen::VectorXd &w) {
  using index_t = Eigen::Index;
  const index_t n = A.cols();
  const index_t m = A.rows();
  Eigen::Map<VectorXd> data( A.data( ), n * m );
  index_t idx = 0;
  assert((n==v.size()) && (n == w.size()));
  for (index_t i=0;i<n;i++) {
    idx += i;
    for (index_t j=i;j<n;j++) {
      w(i) += data[idx]*v(j);
      idx++;
    }
  }
}

void fndriver()
{
  const int n = 10;
  Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> A(n,n);
  Eigen::VectorXd v(n);
  Eigen::VectorXd w(n);
  for (int i=0;i<n;i++) {
    for(int j=0;j<n;j++) {
      A(i,j) = (i+1.0)/(j+1.0);
    }
    v(i) = i+1.0;
  }
  std::cout << "A = " << std::endl << A << std::endl;
  std::cout << "v= [" << v.transpose() << ']' << std::endl;
  w.setZero(); fn(A,v,w);
  std::cout << std::endl << "w = [" << w.transpose() << "]" << std::endl;
}


int main(int argc,char **argv)
{
  const std::span args(argv, argc);
  int code = 0;
  cout << "EIGEN DENSE LINEAR ALGEBRA CODES" << endl;
  if (argc == 2) {
    const int64_t sel = std::strtol(args[1], nullptr, 10);
    switch (sel) {
      case 1: { lsesolve(); break; }
      case 2: { reuselu(); break; }
      case 3: { invpowitdriver(); break; }
      case 4: { fndriver(); break; }
      default: { std::cerr << "Invalid selection" << endl; code = -1; }
    }
  }
  else {
    std::cerr << "Usage: " << args[0] << " <selection>" << endl;
    code = -1;
  }
  return code;
}
