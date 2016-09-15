#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "timer.h"

using namespace Eigen;
using namespace std;

//! \brief Compute the Matrix product $A \times B$ using a naive loop approach (very slow!!)
//! \param[in] A Matrix $2^k \times 2^k$
//! \param[in] B Matrix $2^k \times 2^k$
//! \param[out] Matrix product of A and B of dim $2^k \times 2^k$
MatrixXd loopMatMult(const MatrixXd & A, const MatrixXd & B) 
{
    const int n = A.rows();
    MatrixXd C = MatrixXd::Zero(n,n);

    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          C(i,j) += A(i,k)*B(k,j);
        }
      }
    }

    return C;
}

//! \brief Compute the Matrix product $A \times B$ using Strassen's algorithm.
//! \param[in] A Matrix $2^k \times 2^k$
//! \param[in] B Matrix $2^k \times 2^k$
//! \param[out] Matrix product of A and B of dim $2^k \times 2^k$
MatrixXd strassenMatMult(const MatrixXd & A, const MatrixXd & B)
{
    int n=A.rows();
    MatrixXd C(n,n);
    
    if (n==2)
    {
        C<< A(0,0)*B(0,0) + A(0,1)*B(1,0),
            A(0,0)*B(0,1) + A(0,1)*B(1,1),
            A(1,0)*B(0,0) + A(1,1)*B(1,0),
            A(1,0)*B(0,1) + A(1,1)*B(1,1);
        return C;
    }
    
    else
    {   MatrixXd Q0(n/2,n/2),Q1(n/2,n/2),Q2(n/2,n/2),Q3(n/2,n/2),
        Q4(n/2,n/2),Q5(n/2,n/2),Q6(n/2,n/2);
        
        MatrixXd A11=A.topLeftCorner(n/2,n/2);
        MatrixXd A12=A.topRightCorner(n/2,n/2);
        MatrixXd A21=A.bottomLeftCorner(n/2,n/2);
        MatrixXd A22=A.bottomRightCorner(n/2,n/2);
        
        MatrixXd B11=B.topLeftCorner(n/2,n/2);
        MatrixXd B12=B.topRightCorner(n/2,n/2);
        MatrixXd B21=B.bottomLeftCorner(n/2,n/2);
        MatrixXd B22=B.bottomRightCorner(n/2,n/2);
        
        Q0=strassenMatMult(A11+A22,B11+B22);
        Q1=strassenMatMult(A21+A22,B11);
        Q2=strassenMatMult(A11,B12-B22);
        Q3=strassenMatMult(A22,B21-B11);
        Q4=strassenMatMult(A11+A12,B22);
        Q5=strassenMatMult(A21-A11,B11+B12);
        Q6=strassenMatMult(A12-A22,B21+B22);
        
        C<< Q0+Q3-Q4+Q6 ,
        Q2+Q4,
        Q1+Q3,
        Q0+Q2-Q1+Q5;
        return C;
    }
}


int main(void)
{
    srand((unsigned int) time(0));
    
    //check if loopMatMult and strassenMatMult works
    const int k = 2;
    const int n = pow(2,k);
    MatrixXd A = MatrixXd::Random(n,n),
             B = MatrixXd::Random(n,n),
             AB_loops = loopMatMult(A,B),
             AB_strassen = strassenMatMult(A,B),
             AB_eigen = A*B;
    cout << "Using loop method we get error norm " << (AB_eigen - AB_loops).norm() << endl 
         << "A*B =" << endl << AB_loops << endl;
    cout << "Using Strassen's method we get error norm " << (AB_eigen - AB_strassen).norm() << endl
         << "A*B =" << endl << AB_strassen << endl;
    cout << "Using standard Eigen method A*B =" << endl << AB_eigen << endl;
    
    // compare runtimes of strassenMatMult and loop implementation and Eigens built-in method
    unsigned int repeats = 10;
    std::vector<double> times_eigen, times_strassen, times_loop;
    for(unsigned int k = 4; k <= 10; k++) {
        Timer tm_eigen, tm_strassen, tm_loop;
        for(unsigned int r = 0; r < repeats; ++r) {
            unsigned int n = pow(2,k);
            A = MatrixXd::Random(n,n);
            B = MatrixXd::Random(n,n);
            MatrixXd AB(n,n);
            
            tm_eigen.start();
            AB = A*B;
            tm_eigen.stop();
            
            tm_strassen.start();
            AB = strassenMatMult(A,B);
            tm_strassen.stop();

            tm_loop.start();
            AB = loopMatMult(A,B);
            tm_loop.stop();
        }
        cout << "The loop matrix multiplication took:       " << tm_loop.min() << " s" << endl;
        cout << "The Strassen's algorithm took:       " << tm_strassen.min() << " s" << endl;
        cout << "The standard matrix multiplication took:       " << tm_eigen.min() << " s" << endl;
        
        times_loop.push_back( tm_loop.min() );
        times_strassen.push_back( tm_strassen.min() );
        times_eigen.push_back( tm_eigen.min() );
    }
    
    for(auto it = times_loop.begin(); it != times_loop.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;
    for(auto it = times_eigen.begin(); it != times_eigen.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;
    for(auto it = times_strassen.begin(); it != times_strassen.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;

    return 0;
}
