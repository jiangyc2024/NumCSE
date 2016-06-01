#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <figure.hpp>
#include "../../utils/graphMarketMatrixLoader.hpp"
#include "../../prbuildA/Eigen/prbuildA.hpp"

// sime page rank calculation by tracking fractions of many surfers
void prpowitsim(double d = 0.15, int Nsteps = 5)
{	
	// load connectivity matrix and build transition matrix
	Eigen::MatrixXi G;
	loadGraphMarketMatrix(G, "Harvard500.mtx");
	Eigen::MatrixXd A = prbuildA(G, 0.15);

	// initial distribution
	int N = A.rows();
	Eigen::VectorXd x = Eigen::VectorXd::Ones(N) / N;

	// Plain power iteration for stochastic matrix \Blue{$\VA$}
	for (int l=0; l<Nsteps; ++l)
	{
		x = A*x;
	}

	// plot result
	Eigen::VectorXd pages = Eigen::VectorXd::LinSpaced(N, 1, N);
	mgl::Figure fig;
	std::string title = "bf step " + std::to_string(Nsteps);
	fig.plot(pages, x, " *r");
	fig.title("bf step ");
	fig.xlabel("harvard500: no. of page");
	fig.ylabel("page rank");
	fig.setFontSize(5);
	fig.save("prpowitsim");
}	
