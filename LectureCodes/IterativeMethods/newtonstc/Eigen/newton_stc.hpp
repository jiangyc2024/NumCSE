///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): N.N.
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

/* SAM_LISTING_BEGIN_0 */
template <typename FuncType,typename JacType,typename VecType>
void newton_stc(const FuncType &F,const JacType &DF,
	        VecType &x,double rtol,double atol)
{
  using scalar_t = typename VecType::Scalar;
  scalar_t sn;
  do {
    auto jacfac = DF(x).lu(); // LU-factorize Jacobian \Label[line]{nstc:1}]
    x -= jacfac.solve(F(x));  // Compute next iterate \Label[line]{nstc:2}
    // Compute norm of simplified Newton correction
    sn = jacfac.solve(F(x)).norm();
  }
  // Termination based on simplified Newton correction
  while ((sn > rtol*x.norm()) && (sn > atol));
}
/* SAM_LISTING_END_0 */
