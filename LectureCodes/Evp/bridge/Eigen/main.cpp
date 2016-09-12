#include "trussplot.hpp"
#include <iostream>
#include <Eigen/Dense>


int main()
{
	int n = 11;

	Eigen::MatrixXd pos(n,2);
	Eigen::MatrixXd con(21,2);

	// position of nodes
   	pos << 0,0, 1,0, 2,0, 3,0, 4,0, 5,0, 1,1, 2,1, 3,1, 4,1, 2.5,0.5;

	// indices of connected nodes
	con << 1,2, 2,3, 3,4, 4,5, 5,6, 1,7, 2,7, 3,8, 2,8, 4,8, 5,9, 5,10, 6,10, 7,8, 8,9, 9,10, 8,11, 9,11, 3,11, 4,11, 4,9;
	
	mgl::Figure fig;
	trussplot(fig, pos, con,(char*)"bs-");
	fig.save("bridgetruss");
}
