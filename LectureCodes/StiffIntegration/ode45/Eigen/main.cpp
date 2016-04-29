#include "figure.hpp"
#include "ode45.hpp"
#include <iostream> 


int main()
{
	auto f = [](double x){ return 500*x*x*(1-x); };

	double y0 = 1./100.;
	double T = 1;

	ode45<double> ode = ode45<double>(f);
	ode.options.do_statistics = true;
	ode.options.rtol = 0.001;
	ode.options.atol = 0.0001;

	std::vector<std::pair<double, double>> states = ode.solve(y0, T);
	ode.print(); // prints options and statistics


	// plot states
	//
	
	mgl::Figure lin;

	Eigen::VectorXd y(states.size());
	Eigen::VectorXd t(states.size());
	for (size_t i=0; i<states.size(); ++i)
	{
		y(i) = states[i].first;
		t(i) = states[i].second;
	}

	lin.plot(t, y, "#r^");
	lin.addlabel("y' = 500*y^2*(1-y)", "r");
	lin.legend(1,1);
	lin.title("ODE45");
	lin.xlabel("x");
	lin.ylabel("y");
	lin.setFontSize(3);
	lin.save("ode45");
}
