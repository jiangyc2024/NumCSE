#include "secant.hpp"
#include <cmath>
#include <iostream>


int main()
{
	auto F = [](const double x){ return std::sin(x)*std::cos(x); };
	std::cout << "sin(x)*cos(x) has a zero at " << secant::secant(1, 2, F,1e-6, 1e-6, 100) << std::endl;
}
