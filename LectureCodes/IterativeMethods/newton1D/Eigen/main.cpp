#include <iostream>
#include <cmath>
#include "newton1D.hpp"

int main()
{
	auto F = [](const double x){ return sin(x)*cos(x); };
	auto DF = [](const double x){ return cos(2*x); };
	std::cout << "sin(x)*cos(x) has a zero at " << newton1D(F, DF, 2., 1e-6, 1e-6) << std::endl;
}
