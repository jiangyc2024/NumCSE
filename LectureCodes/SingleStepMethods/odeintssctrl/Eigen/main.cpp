#include "figure.hpp"
#include "odeintssctrl.hpp"
#include <iostream>
#include <math.h>


int main()
{
	auto f = [](double x){ return pow(x,2); };
	auto norm = [](double x){ return fabs(x); };


	// explicit euler (order 1) 
	auto psilow = [&](double h, double y){
		return  y + h*f(y);
	};

	// explicit trapezoidal (order 2)
	auto psihigh = [&](double h, double y){
		double k1 = f(y);
		double k2 = f(y + h*k1);
		return y + (h/2.)*(k1+k2);
	};


	mgl::Figure lin;

	double y0 = 0.5;
	std::vector<std::pair<double, double>> states = odeintssctrl(psilow, 1, psihigh, norm, y0, 1.9, 0.2, 1e-3, 1e-3, 1e-4);

	Eigen::VectorXd y(states.size());
	Eigen::VectorXd t(states.size());
	for (int i=0; i<states.size(); ++i)
	{
		t(i) = states[i].first;
		y(i) = states[i].second;
	}

	lin.plot(t, y, "#r^");
	lin.fplot("1/(2-x)").label("1/(2-x)").style("b-");

	std::cout << "Last state at t="<< states.back().first << " with y=" << states.back().second << std::endl;

	lin.addlabel("y' = y^2", "r");
	lin.legend(1,1);
	lin.title("Refined local stepsize control");
	lin.save("odeintssctrl");


	
}
