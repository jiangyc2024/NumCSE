#include "figure.hpp"
#include "chemstiff.hpp"

void plot(std::vector<std::pair<Eigen::Vector4d, double>> states)
{
	// plot the concentrations of the 4 substances
	//
	
	mgl::Figure lin;

	std::vector<std::vector<double>> y(4, std::vector<double>(states.size()));
	std::vector<double> t(states.size());

	for (size_t i=0; i<states.size(); ++i)
	{
		for (int j=0; j<4; ++j)
		{
			y[j][i] = states[i].first(j);
		}

		t[i] = states[i].second;
	}

	lin.plot(t, y[0]).label("C[a]");
	lin.plot(t, y[1]).label("C[b]");
	lin.plot(t, y[2]).label("C[c]");
	lin.plot(t, y[3]).label("C[d]");

	lin.legend(1,1);
	lin.xlabel("x");
	lin.ylabel("y");
	lin.title("Simulation of 'stiff' chemical reaction");
	lin.setFontSize(3);
	lin.save("chemstiff");
}

int main()
{
	auto states = chemstiff();
	plot(states);	
}
