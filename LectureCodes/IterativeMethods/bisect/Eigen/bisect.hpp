// Searching zero of \Blue{$F$} in \Blue{$[a,b]$} by bisection
template <typename Func>
double bisect(const Func&& F, double a, double b, double tol)
{
	if (a > b)
		std::swap(a,b);

	double fa = F(a), fb = F(b);

	if (fa*fb > 0) 
		throw "f(a) and f(b) have same sign";
	
	int v=1;
	if (fa > 0) v=-1;

	double x = 0.5*(b+a); // determine midpoint
	while (b-a > tol && ((a<x) && (x<b)))
	{
		if (v*F(x) > 0) b=x; // sign(f(x)) == sign(f(b)) -> use x as next right border 
		else a=x; // sign(f(x)) == sign(f(a)) -> use x as next left border 
		x = 0.5*(a+b); // determine next midpoint
	}

	return x;
}
