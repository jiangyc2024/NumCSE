///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <figure/figure.hpp>

#include "spdiags.hpp"

using namespace std;
using namespace Eigen;

int main () {
/* SAM_LISTING_BEGIN_0 */
// Build matrix
int n = 100;
RowVectorXd diag_el(5);	diag_el << -1,-1, 3, -1,-1;
VectorXi diag_no(5);	diag_no << -n, -1, 0, 1, n;
MatrixXd B = diag_el.replicate(2*n,1);
B(n-1,1) = 0; B(n,3) = 0; // delete elements
// A custom function from the Utils folder
SparseMatrix<double> A = spdiags(B, diag_no, 2*n, 2*n);	
// It is not possible to access the LU-factors in the case of
// \eigen's LU-decomposition for sparse matrices.
// Therefore we have to resort to the dense version.
auto solver = MatrixXd(A).lu();
MatrixXd L = MatrixXd::Identity(2*n,2*n);
L += solver.matrixLU().triangularView<StrictlyLower>();
MatrixXd U = solver.matrixLU().triangularView<Upper>();
// Plotting
mgl::Figure fig1, fig2, fig3;
fig1.spy(A);	fig1.setFontSize(4);
fig1.title("Sparse matrix");	fig1.save("sparseA_cpp");
fig2.spy(L);	fig2.setFontSize(4);
fig2.title("Sparse matrix: L factor");	fig2.save("sparseL_cpp");
fig3.spy(U);	fig3.setFontSize(4);
fig3.title("Sparse matrix: U factor");	fig3.save("sparseU_cpp");
/* SAM_LISTING_END_0 */

	return 0;
}
