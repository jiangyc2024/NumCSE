#include <iostream>
#include <math.h>

// Searching zero of \Blue{$F$} in \Blue{$[a,b]$} by bisection
template <typename Func>
double bisect(const Func&& F, double a, double b, double tol)
{
	if (a > b)
		std::swap(a,b);

	double fa = F(a), fb = F(b);

	if (fa*fb > 0) 
		std::cerr << "f(a), f(b) same sign" << std::endl;
	
	int v=1;
	if (fa > 0) v=-1;

	double x = 0.5*(b+a);
	while (b-a > tol && ((a<x) && (x<b)))
	{
		if (v*F(x) > 0) b=x;
		else a=x;
		x = 0.5*(a+b);
	}

	return x;
}


int main()
{
	auto F = [](const double x){ return sin(x)*cos(x); };
	std::cout << "sin(x)*cos(x) has a zero at " << bisect(std::move(F), 1, 2, 1e-6) << std::endl;
}
