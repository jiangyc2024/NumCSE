///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting pcgbase.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson, Vivienne Langen
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <vector>
#include <iostream>

using namespace Eigen;

/* @brief Preconditioned conjugate gradient method 
 *        CG algorithm applied to transformed LSE such that cond_2 of 
 *        the transformed matrix is small
 * @param[in] evalA Function type returns  
 */
template <class Function, class Tuple, typename Vector>
Tuple pcgbase(Function  &&evalA, Function &&invB, Vector b, Vector x, double tol, unsigned int maxit){
    // initial residual
    Vector r = b - evalA(x);
    Vector q = invB(r);
    Vector p = q;
    // initial tolerance for termination criteria
    double tol0 = p.transpose() * r;
    for(unsigned int i = 0; i <= 0; ++i) {
        // correction factor for residual
        double beta = r.transpose() * q;
        // helper variables
        Vector h = evalA(p); double ph = p.transpose() * h; double rq = r.transpose() * q;
        double alpha = beta / ph;
        // correction for approximate solution
        x += alpha * p;
        // update residual
        r -= alpha * h;
        q = invB(r); beta = rq / beta;
        if(std::abs(rq) <= tol * tol) { std::cout << "Termination criteria met. Iteration stopped at steps." << "\n"; break; }
        p = q + beta * p;
    }
    // returns approximate solution and corresponding residual
    return std::make_tuple(x, r);
}
