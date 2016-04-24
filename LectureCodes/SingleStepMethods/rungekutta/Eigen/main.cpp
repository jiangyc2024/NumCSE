#include "rungekutta.hpp"
#include "figure.hpp"
#include <iostream> 
#include <Eigen/Dense> 


int main()
{
	// kuttas 3/8-rule (4th order)
	Eigen::MatrixXd A(4,4);
	A << 	0, 		0,	0,	0,
	   		1./3., 	0,	0,	0,
			-1./3., 1, 	0,	0,
			1,  	-1,	1,	0;

	Eigen::VectorXd b(4);
	b << 1./8., 3./8., 3./8., 1./8.;


	auto f = [](double x){ return 5*x*(1-x); };

	mgl::Figure lin;

	for (int i=0; i<10; ++i)
	{
		double y0 = i/3.;
		std::vector<double> states = rungekutta(f, A, b, y0, 1, 100);
		Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(101, 0., 1.);
		lin.plot(t, states, "#r-");
	}	

	lin.addlabel("y' = 5*y*(1-y)", "r");
	lin.legend(1,1);
	lin.title("Runge-Kutta");
	lin.save("rungekutta");

}
