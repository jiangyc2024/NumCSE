/* SAM_LISTING_BEGIN_1 */
#include <iostream>

#include "strassen.cpp"

int main()
{
#if SOLUTION
    // seed random number generator
    srand((unsigned int) time(0));
    
    // check algorithm for correctness
    int k=2;
    int n=pow(2,k);
    MatrixXd A=MatrixXd::Random(n,n);
    MatrixXd B=MatrixXd::Random(n,n);
    MatrixXd AB(n,n), AxB(n,n);
    AB=strassenMatMult(A,B);
    AxB=A*B;
    std::cout << "Using Strassen's method, A*B=" << AB << std::endl;
    std::cout << "Using standard method, A*B=" << AxB << std::endl;
    std::cout << "The norm of the error is " << (AB-AxB).norm() << std::endl;
#endif
}
/* SAM_LISTING_END_1 */
