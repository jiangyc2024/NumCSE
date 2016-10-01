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
    
    double residual = 0.; // TODO: problem 2d: compute residual
    
    std::cout << "Residual = " << residual << std::endl;
}

