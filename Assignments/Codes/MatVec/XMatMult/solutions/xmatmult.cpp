#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

#include "timer.h"

#include <Eigen/Dense>

using namespace Eigen;

/* @brief Efficiently compute the Matrix-Vector-Product
 * $x=Ay$, where A(i,j) is zero except for values on the
 * diagonal and antidiagonal. There the values are given
 * by the entries of the vector a.
 *
 * @param[in]  a, vector storing the elements of $A$
 * @param[in]  y, vector to multiply $A$ with
 * @param[out] x, result vector of $Ay=x$
 */
 /* SAM_LISTING_BEGIN_1 */
template<class Vector>
void xmatmult(const Vector& a, const Vector& y, Vector& x) {
    assert(a.size() == y.size() && a.size() == x.size()
            && "Input vector dimensions must match");
    unsigned n = a.size();
    
    // looping over half of a,
    // if n is odd we treat the middle element differently
    for (unsigned i = 0; i < n / 2; ++i) {   //integer division
        x(i) = a(i) * y(i) + a(n-i-1) * y(n-i-1);
        x(n-i-1) = x(i);
    }
    
    if (n%2) { //n is odd
        x(n/2) = a(n/2) * y(n/2);
    }
}
/* SAM_LISTING_END_1 */

/*@brief Compares the runtimes of the efficient
 * multipication to the normal Matrix-Vector product
 */
/* SAM_LISTING_BEGIN_2 */
void compare_times() {
    std::cout << "Measuring runtimes for comparison" << std::endl;
    unsigned repeats = 3;
    Timer t_fast, t_slow;
    MatrixXd results(10,3);
    
    for (unsigned k = 14; k > 4; k--) {
        //bitshift operator '<<': 1<<3 == pow(2,3)
        unsigned n = 1<<k;
        VectorXd a,y,x(n);
        a = y = VectorXd::Random(n);
        
        MatrixXd A = a.asDiagonal();
        for (unsigned i = 0; i < n; ++i) {
            A(n-i-1,i) = A(i,i);
        }
        
        //measure multiple times
        for (unsigned i = 0; i < repeats; i++) {
            t_fast.start();
            xmatmult(a,y,x);
            t_fast.stop();
            
            t_slow.start();
            x = A * y;
            t_slow.stop();
        }
        results(k-5,0) = n;
        results(k-5,1) = t_slow.min();
        results(k-5,2) = t_fast.min();
    }
    
    // print results
    std::cout << std::setw(8) << "n"
              << std::setw(15) << "original"
              << std::setw(15) << "efficient"
              << std::setprecision(5) << std::endl;
    for (int i = 0; i < results.rows(); i++) {
        std::cout << std::setw(8)<< results(i,0)
                  << std::setw(15) << results(i,1) << " s"
                  << std::setw(15) << results(i,2) << " s"
                  << std::endl;
    }
}
/* SAM_LISTING_END_2 */

void test() {
    // testing for even n
    unsigned n = 10;
    VectorXd a,y,x(n);
    a = y = VectorXd::Random(n,1);
    //building A for normal Matrix-Vector multiplication O(n*n)
    MatrixXd A = a.asDiagonal();
    for (unsigned i = 0; i < n; ++i) {
        A(n-i-1,i) = A(i,i);
    }
    xmatmult(a,y,x);
    
    double error = (x - A*y).norm();
    if (std::abs(error) > 1e-10) {    //== on floating point numbers is meaninless
        std::cout << "Wrong result for even n" << std::endl;
        return;
    }
    
    // testing for odd n
    n = 11;
    a = y = x= VectorXd::Random(n,1);
    A = a.asDiagonal();
    for (unsigned i = 0; i < n; ++i) {
        A(n-i-1,i) = A(i,i);
    }
    
    xmatmult(a,y,x);
    
    error = (x - A*y).norm();
    if (std::abs(error) > 1e-10) {
        std::cout << "Wrong result for odd n" << std::endl;
        return;
    }
    
    std::cout << "No errors found in your implementation" << std::endl;
}

int main() {
    test();
    compare_times();
}
