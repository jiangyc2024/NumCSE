///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting algcg.m to C++/Eigen.
/// (C) 2021 SAM, D-MATH
/// Author(s): William Andersson
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <iostream>

/* @brief Solves LSE for symmetric positive definite (s.p.d.) matrix applying
 *        the steep linear search onto nested subspaces (Krylov spaces)
 * @param[in] evalA Returns expression $A\times x$ of type as $x$, must pass a
 * handle to a function realising A*x
 * @param[in] b r.h.s. of LSE
 * @param[in] x Initial guess
 * @param[in] tol Tolerance for termination critearium
 * @param[in] maxit Maximum interation step
 * @param[out] x Approxmated solution
 */
template <class Function, class Vector>
Vector cg(Function &&evalA, Vector b, Vector x, double tol,
          unsigned int maxit) {
  // inital residual
  Vector r = b - evalA(x);
  // double rho = 1;
  double n0 = r.norm();
  // first basis vector of Krylov space
  Vector p = r;
  // iterarion to find approx. solution
  for (unsigned int i = 0; i <= maxit; i++) {
    // helper variable
    double rr = r.transpose() * r;
    double beta = rr;
    Vector h = evalA(p);
    double ph = p.transpose() * h;
    double alpha = beta / ph;
    // next guess value
    x += alpha * p;
    // modify residual
    r -= alpha * h;
    // termination of loop approximation sufficient
    if (r.norm() <= tol * n0) {
      std::cout << "Iteration stopped. Termination criteria met. "
                << "\n";
      break;
    }  // relative tolerance (see see Rem.~\ref{rem:cgterm})
    beta = rr / beta;
    p += beta * p;
  }
  return x;
}
