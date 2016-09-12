# include <Eigen/Dense>
using namespace Eigen;

// The function arrowmatvec computes the product the desired product directly as A*A*x
template <class Vector>
void arrowmatvec(const Vector& d, const Vector& a, const Vector& x, Vector& y){
    // Here, we choose a MATLAB style implementation using block construction, you can also use loops
    // If you are interested you can compare both implementation and see if and how they differ
    int n=d.size();
    VectorXd dcut= d.head(n-1);
    VectorXd acut = a.head(n-1);
    MatrixXd ddiag=dcut.asDiagonal();
    MatrixXd A(n,n);
    MatrixXd D = dcut.asDiagonal();
    // If you do not create the temporary matrix D, you will get an error: D must be casted to MatrixXd
    A << D, a.head(n-1), acut.transpose(), d(n-1);
    
    y=A*A*x;
}
