#include <iostream>
#include <cmath>
#include "secant.hpp"

int main()
{
	auto F = [](const double x){ return sin(x)*cos(x); };
	std::cout << "sin(x)*cos(x) has a zero at " << secant(1, 2, std::move(F),1e-6, 1e-6, 100) << std::endl;
}
