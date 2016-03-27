#include <vector>
#include <iostream>
#include <math.h>

template <class Func, class Func2, class NormFunc, class State>
std::vector<std::pair<double, State>> odeintssctrl(Func &Psilow, unsigned int p, Func2 &Psihigh, NormFunc &norm, State& y0, double T, double h0, double reltol, double abstol, double hmin)
{
	double t = 0; // \label{ssctrl:1}
	State y = y0;
	double h = h0;

	std::vector<std::pair<double,State>> states;
	states.push_back({t, y});

	while ((states.back().first < T) && (h >= hmin)) // \label{ssctrl:2}
	{
		State yh = Psihigh(h, y0); // high order discrete evolution \Blue{$\widetilde{\Psibf}^h$}\label{ssctrl:3}
		State yH = Psilow(h, y0); // low order discrete evolution \Blue{${\Psibf}^h$}\label{ssctrl:4}
		double est = norm(yH-yh); // $\leftrightarrow$ \Blue{$\mathrm{EST}_k$}\label{ssctrl:5}
		double tol = std::max(reltol*norm(y0), abstol); // \label{ssctrl:6a}
  
		h = h*std::max(0.5, std::min(2., pow(tol/est, 1./(p+1)))); // Optimal stepsize according to \eqref{eq:ssc}\label{ssctrl:6b}

		if (est < tol) // \label{ssctrl:6}
		{
			//step \Magenta{accepted}\label{ssctrl:7}
			y0 = yh; 
			states.push_back({states.back().first +	std::min(T-states.back().first, h), y0});
		}
	}
	if (h < hmin)
	{
		std::cerr 	<< "Warning: Failure at t=" 
				 	<< states.back().first 
					<< ". Unable to meet integration tolerances without reducing the step size below the smallest value allowed ("<< hmin <<") at time t." << std::endl;
	}

	return states;
}
