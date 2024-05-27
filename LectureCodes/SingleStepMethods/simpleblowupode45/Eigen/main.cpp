#include "ode45.hpp"
#include <iostream> 

//NOLINTBEGIN(bugprone-exception-escape)
int main()
{
	auto f = [](double x){ return pow(x,2); };
	const double y0 = 1;
	const double T = 2;


	Ode45<double> ode = Ode45<double>(f);

	const std::vector<std::pair<double, double>> states = ode.solve(y0, T);

	for (auto state : states)
	{
		std::cout << "t = " << state.first << ", y = " << state.second << std::endl;
	}
}
//NOLINTEND(bugprone-exception-escape)
