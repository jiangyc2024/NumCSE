template <typename FuncType,typename DervType,typename Scalar>
Scalar newton1D(const FuncType &F,const DervType &DF,
		const Scalar &x0,double rtol,double atol)
{
  Scalar s,z = x0;
  do {
    s = F(z)/DF(z);  // compute Newton correction
    z -= s;           // compute next iterate
  }
  // correction based termination (relative and absolute)
  while ((std::abs(s) > rtol*std::abs(z)) && (std::abs(s) > atol));
  return (z);
}
