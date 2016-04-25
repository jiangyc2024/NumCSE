#include "ode45.hpp"
#include <iostream> 


int main()
{
	auto f = [](double x){ return 5*x*(1-x); };
	auto normFunc = [](double x){ return fabs(x); };

	double y0 = 0.2;
	std::vector<std::pair<double, double>> states = ode45(f, y0, 1, normFunc);

	for (auto state : states)
	{
		std::cout << "t = " << state.first << ", y = " << state.second << std::endl;
	}
}
