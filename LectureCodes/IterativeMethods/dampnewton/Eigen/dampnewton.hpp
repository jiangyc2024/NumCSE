#include <cmath>
#include <cstdio>
#include <stdexcept>

template <typename FuncType, typename JacType, typename VecType>
void dampnewton(const FuncType &F, const JacType &DF, VecType &x, double rtol,
                double atol) {
  using index_t = typename VecType::Index;
  using scalar_t = typename VecType::Scalar;
  const index_t n = x.size();
  const scalar_t lmin = 1E-3; // Minimal damping factor
  scalar_t lambda = 1.0;      // Initial and actual damping factor
  VecType s(n);               // Newton corrections
  VecType st(n);
  VecType xn(n);              // Tentative new iterate
  scalar_t sn;                // Norms of Newton corrections
  scalar_t stn;

  do {
    auto jacfac = DF(x).lu(); // LU-factorize Jacobian
    s = jacfac.solve(F(x));   // Newton correction
    sn = s.norm();            // Norm of Newton correction
    lambda *= 2.0;
    do {
      lambda /= 2;
      if (lambda < lmin) {
        throw std::logic_error("No convergence: lambda -> 0");
      }
      xn = x - lambda * s;      // {\bf Tentative next iterate}
      st = jacfac.solve(F(xn)); // Simplified Newton correction
      stn = st.norm();
    } while (stn > (1 - lambda / 2) * sn); // {\bf Natural monotonicity test}
    x = xn;                                // Now: xn accepted as new iterate
    lambda = std::min(2.0 * lambda, 1.0);  // Try to mitigate damping
  }
  // Termination based on simplified Newton correction
  while ((stn > rtol * x.norm()) && (stn > atol));
}
