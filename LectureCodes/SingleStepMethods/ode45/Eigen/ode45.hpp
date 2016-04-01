#include <vector>
#include <Eigen/Dense>
#include <iostream>

template <class Func, class State>
std::vector<State> ode45(Func &f, State &initialState, double T, unsigned int N)
{
	std::vector<State> steps; // saves the state of all steps
	steps.reserve(N+1);
	steps.push_back(initialState);

	double h = T/N;

	Eigen::VectorXd b(7);
	b << 0, 1./5., 3./10., 4./5., 8./9., 1, 1;

	Eigen::MatrixXd A(7,7);
	A << 0,				0,				0,				0,				0,				0,		0,
		 1./5.,			0,				0,				0,				0,				0,		0,
	  	 3./40.,		9./40., 		0,				0,				0,				0,		0,
		 44./45.,		-56./15., 		32./9.,			0,				0,				0,		0,
	 	 19372./6561.,	-25360./2187.,	64448./6561., 	-212./729.,		0,				0,		0,
		 9017./3168.,	-355/33,		46732./5247., 	49./176.,		-5103./18656.,	0,		0,
		 35./384., 		0, 				500./1113.,		125./192.,		-2187./6784.,	11./84.,0;
	
	State y0, y1;
	y0 = initialState;

	for (unsigned int s=0; s < N; ++s)
	{
		y1 = y0;
		std::vector<State> K;
		K.reserve(7);

		for (unsigned int i=0; i<7; ++i)
		{
			// \Vy_{0}+\tstep \sum\limits_{j=1}^{i-1}a_{ij}\Vk_{j})\;
			State incr = y0;
			for (unsigned int j=0; j<i; ++j)
			{
				incr += h*A(i,j)*K.at(j);
			}
			K.push_back(f(incr));

			// \Vy_{1}  += \tstep b_{i}\Vk_{i}\;
			y1 += h*b(i)*K.back();

			steps.push_back(y1);
		}

		y0 = y1;
	}
	return steps;

}
