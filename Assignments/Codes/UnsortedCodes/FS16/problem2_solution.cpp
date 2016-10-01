#include <iostream>

#include <Eigen/Dense>

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

//! \brief Use symmetric Gauss-Seidel iterations to solve the system Ax = b
//! \param[in] A system matrix to be decompsed (L + D + U)
//! \param[in] b r.h.s. vector
//! \param[in,out] x initial guess and last iterate (approximated solution)
//! \param[in] rtol relative tolerance for termination criteria
void GSIt(const Matrix & A, const Vector & b, Vector & x, double rtol) {
    
    // TODO: problem 2c: implement a Gauss-Seidel solver for Ax = b
    
    auto U = Eigen::TriangularView<Matrix, Eigen::StrictlyUpper>(A);
    auto L = Eigen::TriangularView<Matrix, Eigen::StrictlyLower>(A);
    
    auto UpD = Eigen::TriangularView<Matrix, Eigen::Upper>(A);
    auto LpD = Eigen::TriangularView<Matrix, Eigen::Lower>(A);
    
    Vector temp(x.size());
    Vector* xold = &x;
    Vector* xnew = &temp;
    
    unsigned int k = 0;
    double err;
    do {
        *xnew = UpD.solve(b) - UpD.solve(L*LpD.solve(b - U**xold));
        err = (*xold - *xnew).norm();
//         std::cout << k++ << "\t& " << (A*(*xnew) - b).norm() << "\t\\\\" << std::endl;
        std::swap(xold, xnew);
    } while( err > rtol*(*xnew).norm() );
    
    x = *xnew;
    
    return;
}

int main(int, char**) {
    unsigned int n = 9;
    
    Matrix A(n,n);
    for(unsigned int i = 0; i < n; ++i) {
        if(i > 0) A(i,i-1) = 2;
        A(i,i) = 3;
        if(i < n-1) A(i,i+1) = 1;
    }
    Vector b = Vector::Constant(n,1);
    
    //// PROBLEM 2d
    std::cout << "*** PROBLEM 2d:" << std::endl;
    
    // TODO: problem 2d: test the code GSIt using the given data
    Vector x = b;
    GSIt(A, b, x, 10e-8);
    
    double residual = (A*x - b).norm();
    
    std::cout << "Residual = " << residual << std::endl;
}

