#include "../../ode45/Eigen/ode45.hpp"
#include <iostream> 


int main()
{
	auto normFunc = [](double x){ return fabs(x); };

	auto f = [](double x){ return pow(x,2); };
	double y0 = 1;

	std::vector<std::pair<double, double>> states = ode45(f, y0, 2, normFunc);

	for (auto state : states)
	{
		std::cout << "t = " << state.first << ", y = " << state.second << std::endl;
	}
}
