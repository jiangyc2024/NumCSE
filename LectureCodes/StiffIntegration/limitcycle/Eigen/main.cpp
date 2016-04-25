#include "figure.hpp"
#include "limitcycle.hpp"

void plot(mgl::Figure& lin, std::vector<std::pair<Eigen::Vector2d, double>> states)
{
	// plot the solution trajectory
	//

	std::vector<std::vector<double>> y(4, std::vector<double>(states.size()));
	std::vector<double> t(states.size());

	for (size_t i=0; i<states.size(); ++i)
	{
		for (int j=0; j<2; ++j)
		{
			y[j][i] = states[i].first(j);
		}

		t[i] = states[i].second;
	}

	lin.plot(y[0], y[1], "r-");
}

int main()
{
	mgl::Figure lin;
	const double lamda = 5;

	for (int phi=0; phi<360; phi += 20)
	{
		for (int i=0; i<6; ++i)
		{
			double r = i/3.;

			Eigen::Vector2d y0;
			y0 << r*cos(phi), r*sin(phi); // polar to cartesian coordinates

			auto states = limitcycle(lamda, y0);
			plot(lin, states);
		}
	}

	lin.legend(1,1);
	lin.xlabel("y1");
	lin.ylabel("y2");
	lin.title("Strongly attractive limit cycle");
	lin.setFontSize(3);
	lin.save("limitcycle");
}
