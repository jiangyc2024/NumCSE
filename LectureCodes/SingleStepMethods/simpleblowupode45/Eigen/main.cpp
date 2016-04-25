#include "ode45.hpp"
#include <iostream> 


int main()
{
	auto f = [](double x){ return pow(x,2); };
	double y0 = 1;
	double T = 2;


	ode45<double> ode = ode45<double>(f);

	std::vector<std::pair<double, double>> states = ode.solve(y0, T);

	for (auto state : states)
	{
		std::cout << "t = " << state.first << ", y = " << state.second << std::endl;
	}
}
