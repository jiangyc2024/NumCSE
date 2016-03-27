#include <vector>
#include <Eigen/Dense>
#include <iostream>

template <class Func, class State>
std::vector<State> rungekutta(Func &f, const Eigen::MatrixXd &A, const Eigen::VectorXd &b, State &initialState, double T, unsigned int N)
{
	std::vector<State> steps; // saves the state of all steps
	steps.reserve(N+1);
	steps.push_back(initialState);

	double h = T/N;
	
	State y0, y1;
	State *yOld, *yNew;

	// use pointers to avoid copying after each step
	yOld = &y0;
	yNew = &y1;

	y0 = initialState;

	for (unsigned int s=0; s < N; ++s)
	{
		rungekutta_step(f, h, *yOld, *yNew, A, b);
		steps.push_back(y1);
		std::swap(yOld,yNew);
	}

	return steps;
}

template <class Func, class State>
void rungekutta_step(Func f, double h, const State &y0, State &y1, const Eigen::MatrixXd &A, const Eigen::VectorXd &b)
{
	y1 = y0;
	unsigned int s = b.size();

	std::vector<State> K;
	K.reserve(s);

	for (unsigned int i=0; i<s; ++i)
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
	}
}
