#include <Eigen/Dense>

using namespace Eigen;

// The auxiliary function Atimesx computes the function A*x in a smart way, using the particular structure of the matrix A.
template <class Vector>
void Atimesx(const Vector& d, const Vector& a, const Vector& x, Vector& Ax) {
    int n=d.size();
    Ax=(d.array()*x.array()).matrix();
    VectorXd Axcut=Ax.head(n-1);
    VectorXd acut = a.head(n-1);
    VectorXd xcut = x.head(n-1);
    
    Ax << Axcut + x(n-1)*acut, Ax(n-1)+ acut.transpose()*xcut;
}

// We compute A*A*x by using the function Atimesx twice 
template <class Vector>
void arrowmatvec2(const Vector& d, const Vector& a, const Vector& x, Vector& AAx) {
  Eigen::VectorXd Ax;
  Atimesx(d, a, x, Ax);
  Atimesx(d, a, Ax, AAx);
}
