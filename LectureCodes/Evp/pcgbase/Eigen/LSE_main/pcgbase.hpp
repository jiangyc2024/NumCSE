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

using namespace Eigen;

/* @brief Preconditioned conjugate gradient method 
 *        CG algorithm applied to transformed LSE such that cond_2 of 
 *        the transformed matrix is small
 * @param[in] evalA Function, returns A*x, A LVS matrix
 * @param[in] invB Function, returns B*x, B pre-cond. matrix
 * @param[in] b R.h.s.
 * @param[in] x Initial guess
 * @param[in] tol, maxit 
 * @param[out] (x, r) Approx. solution. residual
 */
template <class FunctionA, class FunctionB, typename Vector>
std::pair<Vector, Vector> pcgbase(FunctionA  &&evalA, FunctionB &&invB, Vector b, Vector x, double tol, int maxit){
    // initial residual
    Vector r = b - evalA(x);
    Vector q = invB(r);
    Vector p = q;
    double ry = 1; double ry0 = 1; 
    // this version now is a litterally translation from matlab(instead of pseudo-code)
    for (int i=0; i<maxit; ++i){
        Vector y = invB(r);
        double ry_tmp = ry; 
        ry = r.transpose() * y;
        if(i==1) { p = y; 
                   ry0 = ry; }
        else if(ry < ry0*tol) {
            // if termination criteria met, return solution and residual 
            return std::make_pair(x, r); }
        else {
            double beta = ry / ry_tmp; 
            p = y + beta*p; }
        q = evalA(p); 
        double pq = p.transpose() * q; 
        double alpha = ry / pq;
        // update solution
        x += alpha * p;
        r -= alpha * q;
        } 
    // return solution and residual after maxit
    return std::make_pair(x, r);
}
