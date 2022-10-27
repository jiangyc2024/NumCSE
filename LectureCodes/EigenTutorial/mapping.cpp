// **********************************************************************
// Eigen tutorial codes: http://eigen.tuxfamily.org/dox/GettingStarted.html
// **********************************************************************

#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using std::cout;
using std::endl;

double g(double x) { return x * x; }

// Takes a matrix m and a functor f to apply to all entries of m
template <class Functor> void applyMap(MatrixXd &m, Functor &&f) {
  m = m.unaryExpr(f);
}

void mapping(MatrixXd &m1, const MatrixXd &m2) {
  // Apply a (lambda) function to all entries of a matrix
  const double z(5.0);
  // Capture z
  auto f = [z](double x) { return x + z / x; };
  m1 = m1.unaryExpr(f);

  cout << "f(m1) = " << endl << m1 << endl << endl;

  // Apply an existing, named function (by function pointer) to all entries of a
  // matrix
  m1 = m1.unaryExpr(&g);

  cout << "g(m1) = " << endl << m1 << endl << endl;

  // Apply a binary function to m1 and m2
  auto h = [](double x, double y) { return (x - y) * (x - y); };
  m1 = m1.binaryExpr(m2, h);

  cout << "h(m1, m2) = " << endl << m1 << endl << endl;

  // Pass a lambda function f and m1 to applyMap
  applyMap(m1, f);

  cout << "f(m1) = " << endl << m1 << endl << endl;
}

int main() {
  MatrixXd m1(6, 7);
  MatrixXd m2(6, 7);

  for (int i = 0; i < m1.rows(); i++) {
    for (int j = 0; j < m1.cols(); j++) {
      m1(i, j) = i - j;
      m2(i, j) = j - i;
    }
  }

  cout << "mapping(m1, m2)" << endl << endl;
  mapping(m1, m2);
  return 0;
}
