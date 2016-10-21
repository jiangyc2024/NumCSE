///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Dense>
#include <figure/figure.hpp>

using namespace Eigen;

int main () {
/* SAM_LISTING_BEGIN_0 */
// Study of fill-in with LU-factorization due to pivoting
MatrixXd A(11,11); A.setZero();
A.diagonal() = VectorXd::LinSpaced(11,1,11).cwiseInverse();
A.col(10).setConstant(2); A.row(10).setConstant(2);
auto solver = A.lu();
MatrixXd L = MatrixXd::Identity(11,11);
L += solver.matrixLU().triangularView<StrictlyLower>();
MatrixXd U = solver.matrixLU().triangularView<Upper>();
// Plotting
mgl::Figure fig1, fig2, fig3, fig4;
fig1.spy(A);	fig1.setFontSize(4);
fig1.title("arrow matrix A");	fig1.save("fillinpivotA");
fig2.spy(L);	fig2.setFontSize(4);
fig2.title("L factor");	fig2.save("fillinpivotL");
fig3.spy(U);	fig3.setFontSize(4);
fig3.title("U factor");	fig3.save("fillinpivotU");
std::cout  << A << std::endl;
/* SAM_LISTING_END_0 */
return 0;
}
