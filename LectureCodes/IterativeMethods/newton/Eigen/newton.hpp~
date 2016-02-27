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
