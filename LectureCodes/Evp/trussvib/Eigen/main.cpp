#include "trussvib.hpp"
#include "trussplot.hpp"
#include "eigensolversort.hpp"
#include <iostream>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <string.h>


// plots the truss for the eigenvector with index mode
void plottrussmode(const Eigen::MatrixXd& pos, const Eigen::MatrixXd& con, const Eigen::MatrixXd &eigenvecs, const Eigen::VectorXd &eigenvals, int mode)
{
	int n = pos.rows();
	mode = std::min(n, mode);

	double freq = eigenvals(mode);

	Eigen::MatrixXd u(n,2);
	for (int i=0; i<2*n; ++i)
	{
		if (i%2==0) u(i/2) = eigenvecs(i, mode);
		else u((i-1)/2,1) = eigenvecs(i, mode);
	}

	mgl::Figure figMode;
	trussplot(figMode, pos,con, (char*)"ks--");
	trussplot(figMode, pos + 0.5*u, con, (char*)"bs-");
	figMode.title("mode " + std::to_string(mode) + ": frequency = " + std::to_string(freq));
	figMode.save(std::string("trussmode") + std::to_string(mode));
}

int main()
{
	int n = 11;

	Eigen::MatrixXd pos(n,2);
	Eigen::MatrixXd con(21,2);

	// position of nodes
   	pos << 0,0, 1,0, 2,0, 3,0, 4,0, 5,0, 1,1, 2,1, 3,1, 4,1, 2.5,0.5;

	// indices of connected nodes
	con << 1,2, 2,3, 3,4, 4,5, 5,6, 1,7, 2,7, 3,8, 2,8, 4,8, 5,9, 5,10, 6,10, 7,8, 8,9, 9,10, 8,11, 9,11, 3,11, 4,11, 4,9;
	

	// Visualizing nonzero eigenmodes of an elastic truss

	// Compute eigenmodes
	Eigen::EigenSolver<Eigen::MatrixXd> ev = trussvib(pos, con);
	Eigen::MatrixXd eigenvecs;
	Eigen::VectorXd eigenvals;
	std::tie(eigenvals, eigenvecs) = eigensolversort(ev);

	// Plot resonances
	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(2*n,1,2*n);
	mgl::Figure fig;
	fig.plot(x,eigenvals, "r+");
	fig.xlabel("no. of eigenvalue}");
	fig.ylabel("eigenvalue");
	fig.title("Truss resonant frequencies");
	fig.save("trussfreq");

	for (int i=3; i<=6; ++i)
	{
		plottrussmode(pos, con, eigenvecs, eigenvals, i);
	}
}
