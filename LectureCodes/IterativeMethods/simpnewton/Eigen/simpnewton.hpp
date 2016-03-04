template<typename Func, typename Jac, typename Vec>
void simpnewton(Vec& x, Func F, Jac DF, double rtol, double atol)
{

	// C++ template for simplified Newton method
	auto lu = DF(x).lu(); // LU decomposition
	Vec s;
	double ns,nx;

	do
	{
		s = lu.solve(F(x));
		x = x-s;
		ns = s.norm();
		nx = x.norm();
	}
	// termination based on relative and absolute tolerance
	while((ns > rtol*nx) && (ns > atol));
}
