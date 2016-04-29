#include "Eigen/Dense"
#include "../../ode45/Eigen/ode45.hpp"

std::vector<std::pair<Eigen::Vector4d, double>> chemstiff()
{
	// Simulation of kinetics of coupled chemical reactions with vastly different reaction
	// rates, see \eqref{eq:chemstiff} for the ODE model.
	// reaction rates \Blue{$k_1,k_2,k_3,k_4$}, \Magenta{$k_1,k_2 \gg k_3,k_4$}.
	const double k1 = 1E4, k2 = 1E3, k3 = 10, k4 = 1;

	// definition of right hand side function for ODE solver
	auto f = [&](const Eigen::Vector4d y) {

		Eigen::Vector4d res;
		res << 	-k1*y(0)*y(1) + k2*y(2) - k3*y(0)*y(2) +  k4*y(3),
				-k1*y(0)*y(1) + k2*y(2),
				k1*y(0)*y(1) - k2*y(2) - k3*y(0)*y(2) + k4*y(3),
				k3*y(0)*y(2) - k4*y(3);

		return res;
	};

	const double T = 1; // Integration time interval

	Eigen::Vector4d y0; // Initial value \Blue{$\Vy_0$}
	y0 << 1, 1, 10, 0; 

	// compute ``exact'' solution, using \texttt{ode113} with tight error tolerances
	ode45<Eigen::Vector4d> ode = ode45<Eigen::Vector4d>(f);
	ode.options.do_statistics = true;

	auto norm = [](const Eigen::Vector4d v) { return v.norm(); };

	std::vector<std::pair<Eigen::Vector4d, double>> states = ode.solve(y0, T, norm);
	ode.print(); // prints options and statistics

	return states;
}
