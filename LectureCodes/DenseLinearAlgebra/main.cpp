// **********************************************************************
// Eigen codes for testing dense linear algebra routines
// **********************************************************************

#include <assert.h>
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void lsesolve(std::size_t n = 7)
{
  // Initialize a special invertible matrices
  MatrixXd mat = MatrixXd::Identity(n,n) +
    VectorXd::Constant(n,1.0)*RowVectorXd::Constant(n,1.0);
  cout << "Matrix mat = " << endl << mat << endl;
  MatrixXd utm = mat.triangularView<Upper>();
  cout << "Matrix utm = " << endl << utm << endl;
  // Multiple right hand side vector
  MatrixXd B = MatrixXd::Random(n,2);
  // Solve linear system using various decompositions
  MatrixXd X = mat.lu().solve(B);
  MatrixXd X2 = mat.fullPivLu().solve(B);
  MatrixXd X3 = mat.householderQr().solve(B);
  MatrixXd X4 = mat.llt().solve(B);
  MatrixXd X5 = mat.ldlt().solve(B);
  cout << "|X2-X| = " << (X2-X).norm() << endl;
  cout << "|X3-X| = " << (X3-X).norm() << endl;
  cout << "|X4-X| = " << (X4-X).norm() << endl;
  cout << "|X5-X| = " << (X5-X).norm() << endl;
}

void reuselu(double tol = 1E-3,std::size_t n = 7)
{
  // Initialize a special invertible matrices
  MatrixXd A = MatrixXd::Identity(n,n) +
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
  using scalar_t = typename VecType::Scalar;
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

void invpowitdriver(double tol = 1E-3,std::size_t n = 7)
{
  // Initialize a special invertible matrices
  MatrixXd A = MatrixXd::Identity(n,n) +
    VectorXd::Constant(n,1.0)*RowVectorXd::Constant(n,1.0);
  VectorXd ev = invpowit<VectorXd>(A,tol);
  cout << "ev = [" << ev.transpose() << "]" << endl;
  VectorXd av = A*ev;
  cout << "ev/av = [" << (ev.array()/av.array()).transpose() << "]" << endl;
}

void fn(const Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> &A,
	const Eigen::VectorXd &v,Eigen::VectorXd &w) {
  using scalar_t = double;
  using index_t = std::size_t;
  const scalar_t *a = A.data();
  const index_t n = A.rows();
  assert((n==v.size()) && (n == w.size()));
  for (index_t i=0;i<n;i++) {
    a += i;
    for (index_t j=i;j<n;j++) {
      w(i) += (*a)*v(j);
      a++;
    }
  }
}

void fndriver(void)
{
  const int n = 10;
  Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> A(n,n);
  Eigen::VectorXd v(n),w(n);
  for (int i=0;i<n;i++) {
    for(int j=0;j<n;j++) A(i,j) = (i+1.0)/(j+1.0);
    v(i) = i+1.0;
  }
  std::cout << "A = " << std::endl << A << std::endl;
  std::cout << "v= [" << v.transpose() << ']' << std::endl;
  w.setZero(); fn(A,v,w);
  std::cout << std::endl << "w = [" << w.transpose() << "]" << std::endl;
}


int main(int argc,char **argv)
{
  cout << "EIGEN DENSE LINEAR ALGEBRA CODES" << endl;
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <selection>" << endl;
    return(-1L);
  }
  else {
    const int sel = atoi(argv[1]);
    switch (sel) {
    case 1: { lsesolve(); break; }
    case 2: { reuselu(); break; }
    case 3: { invpowitdriver(); break; }
    case 4: { fndriver(); break; }
    default: { cerr << "Invalid selection" << endl; exit(-1L); }
    }
  }
  return(0);
}
