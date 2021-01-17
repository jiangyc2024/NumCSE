///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting gradit.m to C++/Eigen.
/// (C) 2021 SAM, D-MATH
/// Author(s): William Andersson, Vivienne Langen
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

/* @brief Solves LSE for symmetric positive definite (s.p.d.) matrix using
 *        Lemma 10.1.1.1 and the method of steepest descent, see Section 10.1.3
 * @param[in] evalA Returns expression $A\times x$ of type as $x$
 * @param[in] A Quadratic, symmetric, positive definite matrix
 * @param[in] b Vector type, r.h.s. of LSE
 * @param[in] x Vector type, initial guess
 * @param[in] rtol, atol Relative, absolute tolerance for termination criteria
 * @param[in] maxit Maximum iteration steps for method
 * @param[out] x Approximate solution to LSE
 */
template <class Function, typename Matrix, typename Vector>
Vector gradit(Function&& evalA, Matrix A, Vector b, Vector x, double rtol,
              double atol, unsigned int maxit) {
  // inital residual r0
  Vector r_k = b - evalA(x);
  // iteration to find x_approx
  for (unsigned int k = 0; k <= maxit; k++) {
    // variable to minimize evaluation
    Vector Ar_k = evalA(r_k);
    // paramter for correction
    double a1 = r_k.transpose() * r_k;
    double a2 = r_k.transpose() * Ar_k;
    double tstar = a1 / a2;  // scalar
    // update approximate solution
    x += tstar * r_k;
    // norm of correction
    // note that: x_k+1 = x_k + tstar * r_k
    //            => x_k+1 - x_k = tstar * r_k
    // thus valid as termination criteria
    double cn = std::abs(tstar) * r_k.norm();
    if (cn < rtol * x.norm() || cn < atol) {
      std::cout << "Iteration stopped. Termination criteria met."
                << "\n";
      break;
    };
    // updating the residual
    r_k -= tstar * Ar_k;
  }
  return x;
}
