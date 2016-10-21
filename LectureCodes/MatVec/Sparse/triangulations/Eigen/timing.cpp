///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Thomas Etterlin <thomaset@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

#include "triangulation.hpp"
#include "timer.h"

int main(){
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
Eigen::VectorXd x_ref, y_ref, xs, ys;
Eigen::MatrixXi T_ref;
int refine_steps = 8;
int nRuns = 3;
Eigen::MatrixXd result(refine_steps, 6);
typedef Eigen::SparseMatrix<double> SpaMat;
// Timing of solving graph laplacian
for(int i = 1; i <= refine_steps; ++i){
	Timer t;
	t.start();
	refinemesh(x, y, T, x_ref, y_ref, T_ref);
	t.stop();
	Eigen::SparseMatrix<double> A_int;
	Eigen::MatrixXd rhs;
	smoothmesh_analysis(x_ref, y_ref, T_ref, xs, ys, A_int, rhs);
	// Timing
	Timer t1, t2, t3, t4, t5;
	for(int k = 0; k < nRuns; ++k){
		// ############## SparseLU ##############
		t1.start();
		Eigen::SparseLU<SpaMat> solver1;
		// Indicate that the pattern of the input matrix is symmetric
		solver1.isSymmetric(true);
		// Compute decomposition
		solver1.compute(A_int);
		// Solve the linear system of equations with the corresponding rhs
		Eigen::MatrixXd xy_int = solver1.solve(rhs);
		t1.stop();
		// ############## SimplicialLDLT ##############
		t2.start();
		Eigen::SimplicialLDLT<SpaMat> solver2(A_int);
		xy_int = solver2.solve(rhs);
		t2.stop();
		// ############## ConjugateGradient ##############
		t3.start();
		Eigen::ConjugateGradient<SpaMat> solver3(A_int);
		// Use coordinates after refinement as guess instead of 0
		Eigen::MatrixXd guess(x_ref.size(),2);
		guess << x_ref, y_ref;
		xy_int = solver3.solveWithGuess(rhs, guess);
		t3.stop();
		// ############## PardisoLU ##############
		t4.start();
		Eigen::PardisoLU<SpaMat> solver4(A_int);
		xy_int = solver4.solve(rhs);
		t4.stop();
		// ############## PardisoLDLT ##############
		t5.start();
		Eigen::PardisoLDLT<SpaMat> solver5(A_int);
		xy_int = solver5.solve(rhs);
		t5.stop();
	}
	result(i-1, 0) = A_int.rows();
	result(i-1,1) = t1.min();
	result(i-1,2) = t2.min();
	result(i-1,3) = t3.min();
	result(i-1,4) = t4.min();
	result(i-1,5) = t5.min();
	x = xs;
	y = ys;
	T = T_ref;
	Eigen::MatrixXi E, Eb;
	// Extract the edge information of a mesh
	// E and Eb are matrices whose rows contain the numbers of the 
	// endpoints of edges
	processmesh(T,E,Eb);		
	int Nv = T.maxCoeff()+1;
	int Ne = E.rows();
	int Nt = T.rows();
}
std::cout << std::scientific << result << std::endl;
return 0;
}
