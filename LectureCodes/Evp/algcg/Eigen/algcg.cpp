///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting algcg.m to C++/Eigen.
/// (C) 2020 SAM, D-MATH
/// Author(s): William Andersson
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>

template <class Function, class Vector>
Vector cg(Function &&evalA, Vector b, Vector x, double tol,
          unsigned int maxit) {
  // x supplies the initial guess,
  // maxit the maximal number of CG steps.
  // evalA must pass a handle to a function realising A*x

  Vector r = b - evalA(x);
  double rho = 1;
  double n0 = r.norm();
  Vector p;

  for (unsigned int i = 0; i < maxit; i++) {
    double rho_old = rho;
    rho = r.transpose() * r;

    if (i == 0) {
      p = r;
    }

    else {
      double beta = rho / rho_old;
      p = r + beta * p;
    }

    Vector q = evalA(p);
    double alpha = rho / (p.transpose() * q);

    x += alpha * p;  // Update approximate solution

    // Termination, see Rem. cgterm
    if ((b - evalA(x)).norm() <= tol * n0) {
      // ^ MAYBE the A * x should be evalA(x)? Since
      // evalA(x) should "realise A * x"?
      break;
    }

    r -= alpha * q;  // update of residual
  }

  return x;
}
