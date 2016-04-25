#include "Eigen/Dense"
#include "../../ode45/Eigen/ode45.hpp"

std::vector<std::pair<Eigen::Vector2d, double>> limitcycle(double lamda, Eigen::Vector2d y0)
{
	// define right hand side vectorfield
	auto f = [&](const Eigen::Vector2d y) {

		Eigen::Vector2d res;
		res << 	-y(1), y(0);
		res += lamda*(1-y(0)*y(0) - y(1)*y(1))*y;

		return res;
	};

	const double T = 2*M_PI; // Integration time interval

	ode45<Eigen::Vector2d> ode = ode45<Eigen::Vector2d>(f);
	ode.options.do_statistics = true;

	// solve with ode45
	std::vector<std::pair<Eigen::Vector2d, double>> states = ode.solve(y0, T);

	return states;
}
