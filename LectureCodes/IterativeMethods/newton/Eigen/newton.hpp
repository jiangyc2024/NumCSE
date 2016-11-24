///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair (hiptmair@sam.math.ethz.ch)
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

/* SAM_LISTING_BEGIN_0 */
template <typename FuncType,typename JacType,typename VecType>
void newton(const FuncType &F,const JacType &DFinv,
	    VecType &x,double rtol,double atol)
{
  using index_t = typename VecType::Index;
  using scalar_t = typename VecType::Scalar;
  const index_t n = x.size();
  VecType s(n);
  scalar_t sn;
  do {
    s = DFinv(x,F(x));  // compute Newton correction
    x -= s;                  // compute next iterate
    sn = s.norm();
  }
  // correction based termination (relative and absolute)
  while ((sn > rtol*x.norm()) && (sn > atol));
}
/* SAM_LISTING_END_0 */
