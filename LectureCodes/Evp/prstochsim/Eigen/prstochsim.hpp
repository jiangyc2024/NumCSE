#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <figure.hpp>
#include <string>
#include <stdlib.h>
#include <time.h>
#include "../../utils/graphMarketMatrixLoader.hpp"


// simplified page rank algorithm: simulates Nhops link clicks of a surfer
void prstochsim(std::string path, int Nhops)
{

	// Load web graph data stored in \Blue{$N\times N$}-matrix \Blue{$\VG$}
	Eigen::SparseMatrix<int> G;
	if (!loadGraphMarketMatrix<int>(G, path))
	{
		std::cerr << "Matrix Hardvard500.mtx has not been found." << std::endl;
		exit(EXIT_FAILURE);
	}

	int N = G.rows();
	double d = 0.15;
	
	int cp = 0; // current page index
	Eigen::VectorXd count(N); // how often each page was visited

	srand(time(NULL)); // initialize random seed

	for (int i=0; i<Nhops; ++i)
	{
		// Find links from current page \texttt{cp}
		Eigen::VectorXd col = Eigen::VectorXd(G.col(cp));
		std::vector<int> indices;	
		for (int j=0; j<N; ++j) if (col(j) != 0) indices.push_back(j);


		double rn = ((double) rand() / (RAND_MAX)); // random value in [0,1]

	  	// If no links, jump to any other pages with equal probability
		if (indices.size() == 0)	cp = floor(rn*N);

		// With probabilty \Blue{$d$} jump to any other page  
		else if (rn < d) 			cp = floor(rn/d*N);

	  	// Follow outgoing links with equal probabilty
		else cp = indices[floor((rn-d)/(1.-d)*indices.size())];

		count(cp)++;
	}	
	
	// plot result
	Eigen::VectorXd pages = Eigen::VectorXd::LinSpaced(N, 1, N);
	mgl::Figure fig;
	count /= Nhops; // normalize visits
	fig.plot(pages, count, " .r");
	fig.xlabel("harvard500: no. of page");
	fig.ylabel("page rank");
	std::string title = "page rank, harvard500: " + std::to_string(Nhops) + " hops";
	fig.title(title);
	fig.setFontSize(5);
	fig.save("prstochsim");
}
