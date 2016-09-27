#include <iostream>

#include <Eigen/Dense>

#include "timer.h"

using namespace Eigen;

/* @brief Compute lmin from vector d, naive implementation
 * @param[in] d An n-dimensional vector
 * @param[in] tol Scalar of type 'double', the tolerance
 * @param[out] lmin Scalar of type 'double'
 */
/* SAM_LISTING_BEGIN_0 */
void rankoneinvit(const VectorXd & d, const double & tol, double & lmin)
{
    VectorXd ev = d;
    lmin = 0;
    double lnew = d.cwiseAbs().minCoeff();

    while(abs(lnew-lmin)>tol*lmin) {
#if SOLUTION
        tm_slow.start();
        lmin = lnew;
        MatrixXd M = d.asDiagonal();
        M += ev*ev.transpose();
        ev = M.lu().solve(ev);
        ev.normalize();
        lnew = ev.transpose()*M*ev;
        tm_slow.stop();
#endif // SOLUTION
    }

    lmin = lnew;
}
/* SAM_LISTING_END_0 */

/* @brief Compute lmin from vector d, optimized implementation
 * @param[in] d An n-dimensional vector
 * @param[in] tol Scalar of type 'double', the tolerance
 * @param[out] lmin Scalar of type 'double'
 */
/* SAM_LISTING_BEGIN_1 */    
void rankoneinvit_fast(const VectorXd & d, const double & tol, double & lmin)
{
    VectorXd ev=d;
    lmin=0;
    double lnew=d.cwiseAbs().minCoeff();
    
    VectorXd dinv=(1/d.array()).matrix();
    while (abs(lnew-lmin)>tol*lmin) {
#if SOLUTION
        tm_fast.start();

        lmin = lnew;
        VectorXd ev0 = ev;
        
	// Here we solve the linear system
	// with the Sherman-Morrison-Woodbury formula
	// in the case of rank-1 perturbations
        VectorXd Aib = dinv.cwiseProduct(ev);
        double temp = ev.transpose()*Aib;
        ev = Aib*(1-temp/(1+temp));
        
        ev.normalize();
	// Better than the corresponding naive implementation
        lnew = ev.transpose()*d.cwiseProduct(ev) + pow(ev.transpose()*ev0,2);

        tm_fast.stop();
#endif // SOLUTION
    }

    lmin=lnew;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */ 
int main() {
    srand((unsigned int) time(0));
    double tol = 1e-3;
    double lmin;
    int n = 10;
    
    // Compute with both implementations
    VectorXd d = VectorXd::Random(n);
    std::cout << "Direct porting from MATLAB (naive implementation): " << std::endl;
    rankoneinvit(d,tol,lmin);
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "Fast implementation: " << std::endl;
    rankoneinvit_fast(d,tol,lmin);
    std::cout << "lmin = " << lmin << std::endl;
    
    // Compare runtimes of different implementations of rankoneinvit
    std::cout << "*** Runtime comparison of two implementations" << std::endl;
    unsigned int repeats = 3;
    timer<> tm_slow, tm_fast;
    
    for(unsigned int p = 2; p <= 9; p++) {
        tm_slow.reset();
        tm_fast.reset();
        unsigned int n = pow(2,p);
        
        for(unsigned int r = 0; r < repeats; ++r) { 
         // d = VectorXd::Random(n);
            d = VectorXd::LinSpaced(n,1,2);

	    tm_slow.start();
            rankoneinvit(d,tol,lmin);
	    tm_slow.stop();
            
	    tm_fast.start();
            rankoneinvit_fast(d,tol,lmin);
	    tm_fast.stop();
        }
        
        std::cout << "The slow method took: " << tm_slow.min().avg().count() / 1000000. << " ms for n = " << n << std::endl;
        std::cout << "The fast method took: " << tm_fast.min().avg().count() / 1000000. << " ms for n = " << n << std::endl;
    }
}
/* SAM_LISTING_END_2 */
