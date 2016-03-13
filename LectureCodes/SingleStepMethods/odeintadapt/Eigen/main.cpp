#include "odeintadapt.hpp"
#include <iostream>
#include <math.h>


int main()
{
	auto f = [](double x){ return pow(x,2); };
	auto norm = [](double x){ return fabs(x); };


	// explicit euler (order 1) 
	auto psilow = [&](double h, double y){
		return  y + h*f(y);
	};

	// explicit trapezoidal (order 2)
	auto psihigh = [&](double h, double y){
		double k1 = f(y);
		double k2 = f(y + h*k1);
		return y + (h/2.)*(k1+k2);
	};

	double y0 = 0.5;
	std::vector<std::pair<double, double>> states = odeintadapt(psilow, psihigh, norm, y0, 2., 0.1, 1e-6, 1e-6, 1e-8);

	std::cout << "Last state at t="<< states.back().first << " with y=" << states.back().second << std::endl;
	
}
