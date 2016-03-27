#include "ode45.hpp"
#include <iostream> 


int main()
{
	auto f = [](double x){ return 5*x*(1-x); };
	double y0 = 10;
	std::vector<double> states = ode45(f, y0, 1, 100);

	std::cout << states.back() << std::endl;
}
