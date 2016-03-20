#include "../../../SingleStepMethods/ode45/Eigen/ode45.hpp"
#include <iostream> 


int main()
{
	auto f = [](double x){ return 500*x*x*(1-x); };
	auto normFunc = [](double x){ return fabs(x); };

	double y0 = 0.01;
	std::vector<std::pair<double, double>> states = ode45(f, y0, 1, normFunc, 1e-2, 1e-1, 1e-2, 1e-5, false);

	for (auto state : states)
	{
		std::cout << "t = " << state.first << ", y = " << state.second << std::endl;
	}
}
