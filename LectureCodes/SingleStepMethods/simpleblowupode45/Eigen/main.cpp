#include "../../ode45/Eigen/ode45.hpp"
#include <iostream> 


int main()
{
	auto f = [](double x){ return pow(x,2); };
	double y0 = 1;
	std::vector<double> states = ode45(f, y0, 2, 10);

	for (double y : states)
		std::cout << y << std::endl;
}
