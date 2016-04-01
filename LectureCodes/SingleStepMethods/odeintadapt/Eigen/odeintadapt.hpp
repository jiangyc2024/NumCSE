#include <vector>
#include <iostream>

template <class Func, class Func2, class NormFunc, class State>
std::vector<std::pair<double, State>> odeintadapt(Func &Psilow, Func2 &Psihigh, NormFunc &norm, State& y0, double T, double h0, double reltol, double abstol, double hmin)
{
	double t = 0; // \label{odeintadapt:1}
	State y = y0;
	double h = h0;

	std::vector<std::pair<double,State>> states;
	states.push_back({t, y});

	while ((states.back().first < T) && (h >= hmin)) // \label{odeintadapt:2}
	{
		State yh = Psihigh(h, y0); // high order discrete evolution \Blue{$\widetilde{\Psibf}^h$} \label{odeintadapt:3}
		State yH = Psilow(h, y0); // low order discrete evolution \Blue{${\Psibf}^h$} \label{odeintadapt:4}
		double est = norm(yH-yh); // $\leftrightarrow$ \Blue{$\mathrm{EST}_k$}\label{odeintadapt:5}

		if (est < std::max(reltol*norm(y0), abstol)) // \label{odeintadapt:6}
		{
			y0 = yh; 
			states.push_back({states.back().first +	std::min(T-states.back().first, h), y0}); // \label{odeintadapt:7}
			h = 1.1*h;  // step \Magenta{accepted}, try with increased stepsize \label{odeintadapt:8}
		}
		else
		{
			h = h/2; // step \Magenta{rejected}, try with half the stepsize \label{odeintadapt:9}
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
