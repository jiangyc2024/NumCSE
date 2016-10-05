#include <iostream>

#include <Eigen/Dense>

#include <figure/figure.hpp>
#include "triangulation.hpp"



// code is ported from Matlab --> not very inteligent - yet?

int main(){
/* SAM_LISTING_BEGIN_0 */
//	Demonstration for visualizing a plane triangular mesh
// Initialize node coordinates
Eigen::VectorXd x(10), y(10);
// x and y coordinate of mesh
x << 1.0,0.60,0.12,0.81,0.63,0.09,0.27,0.54,0.95,0.96;
y << 0.15,0.97,0.95,0.48,0.80,0.14,0.42,0.91,0.79,0.95;
// Then specify triangles through the indices of their vertices. These
// indices refer to the ordering of the coordinates as given in the
// vectors x and y.
Eigen::MatrixXi T(11,3);
T << 7, 1, 2,   5, 6, 2,    4, 1, 7,    6, 7, 2,
	6, 4, 7,   6, 5, 0,    3, 6, 0,    8, 4, 3, 
	3, 4, 6,   8, 1, 4,    9, 1, 8;
// Call the \figureclass plotting routine, draw mesh with blue edges
// red vertices and a numbering/ordering
mgl::Figure fig1; fig1.setFontSize(8);
fig1.ranges(0.0, 1.05, 0.0, 1.05);
fig1.triplot(T, x, y, "b?"); // drawing triangulation with numbers
fig1.plot(x, y, " *r"); // mark vertices
fig1.save("meshplot_cpp");
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd x_ref, y_ref, xs, ys;
Eigen::MatrixXi T_ref;
	int refine_steps = 3;
	for(int i = 1; i <= refine_steps; ++i){
		refinemesh(x, y, T, x_ref, y_ref, T_ref);
		mgl::Figure fig2;
		fig2.title("Refined mesh level " + std::to_string(i));
		fig2.ranges(0.0, 1.05, 0.0, 1.05);
		fig2.triplot(T_ref, x_ref, y_ref, "b");
		fig2.save("rmesh" + std::to_string(i) + "_cpp");
		smoothmesh(x_ref, y_ref, T_ref, xs, ys);
		x = xs;
		y = ys;
		T = T_ref;
		mgl::Figure fig3;
		fig3.title("Smoothed mesh level " + std::to_string(i));
		fig3.ranges(0.0, 1.05, 0.0, 1.05);
		fig3.triplot(T, x, y, "b");
		fig3.save("smesh" + std::to_string(i) + "_cpp");
	}
/* SAM_LISTING_END_1 */	
	return 0;
}
