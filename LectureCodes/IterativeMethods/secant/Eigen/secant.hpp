
template<typename Func>
double secant(double x0, double x1, Func&& F, double rtol, double atol, unsigned int maxIt)
{
	double fo = F(x0);
	for (unsigned int i=0; i<maxIt; ++i)
	{
		double fn = F(x1);
		double s = fn*(x1-x0)/(fn-fo); // secant correction
		x0 = x1; 
		x1 = x1-s;

		// correction based termination (relative and absolute)
		if (std::abs(s) < std::max(atol,rtol*std::min(std::abs(x0),std::abs(x1))))
			return x1;

		fo = fn; 
	}
	return x1;
}
