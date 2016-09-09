#include <iostream>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra> // for import
#include <figure/figure.hpp>

using namespace Eigen;

int main (int argc, char* argv[]) {
	typedef SparseMatrix<double> SpMat_t;
	SpMat_t M;
	// load file
	std::string filename = argv[1];
	if(loadMarket(M,filename) != true){
		std::cout << "failed import" << std::endl;
	}
/* SAM_LISTING_BEGIN_0 */
// L and U cannot be extracted from SparseLU --> LDLT
SimplicialLDLT<SpMat_t, Lower, AMDOrdering<int> > solver1(M);
SimplicialLDLT<SpMat_t, Lower, NaturalOrdering<int> > solver2(M);
MatrixXd U1 = MatrixXd(solver1.matrixU());
MatrixXd U2 = MatrixXd(solver2.matrixU());
// Plotting
mgl::Figure fig1, fig2, fig3;
fig1.spy(M);	fig1.setFontSize(4); fig1.save("Mspy");
fig2.spy(U1);	fig2.setFontSize(4); fig2.save("AMDUSpy");
fig3.spy(U2);	fig3.setFontSize(4); fig3.save("NaturalUSpy");
/* SAM_LISTING_END_0 */
	return 0;
}
