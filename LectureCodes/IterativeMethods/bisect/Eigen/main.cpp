#include <iostream>
#include <math.h>
#include "bisect.hpp"

int main()
{
	auto F = [](const double x){ return sin(x)*cos(x); };
	std::cout << "sin(x)*cos(x) has a zero at " << bisect(std::move(F), 1, 2, 1e-6) << std::endl;
}
