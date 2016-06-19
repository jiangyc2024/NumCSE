#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <figure/figure.hpp>
using namespace std;
using namespace Eigen;


int main () {
	typedef SparseMatrix<double> SpMat_t;
	SpMat_t M;
	// load file
	std::string filename = "../poisson2D.mtx";
	if(loadMarket(M,filename) != true){
		std::cout << "failed import" << std::endl;
	}
	#pragma begin<0>
	// L and U cannot be extracted from SparseLU --> LDLT
	SimplicialLDLT<SpMat_t, Lower, AMDOrdering<int> > solver1(M);
	SimplicialLDLT<SpMat_t, Lower, NaturalOrdering<int> > solver2(M);
	MatrixXd U1 = MatrixXd(solver1.matrixU());
	MatrixXd U2 = MatrixXd(solver2.matrixU());
	// Plotting
	mgl::Figure fig1, fig2, fig3;
	fig1.spy(M);	fig1.setFontSize(4); fig1.grid(false);
	fig1.title("arrow matrix A");	fig1.save("Mspy");
	fig2.spy(U1);	fig2.setFontSize(4); fig2.grid(false);
	fig2.title("U factor");	fig2.save("AMDMSpy");
	fig3.spy(U2);	fig3.setFontSize(4); fig3.grid(false);
	fig3.title("U factor");	fig3.save("NaturalMSpy");
	#pragma end<0>
	return 0;
}
