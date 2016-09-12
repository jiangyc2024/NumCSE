#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include "spdiags.hpp"
#include "imgsegmat.hpp"
#include "imread.hpp"

// segmentation of P: build matrices, see \cref{mc:imgsegmat} and \eqref{eq:segmat}
void seg(Eigen::MatrixXd &P)
{
	int m = P.rows(); int n = P.cols();	

	Eigen::SparseMatrix<double> D;
	Eigen::SparseMatrix<double> A;

	auto Sfun = [](double x, double y) 
		{ return std::exp(-0.1*std::pow(x-y, 2)); };

	std::tie(A,D) = imgsegmat(P, Sfun);

	// Build scaling matrics 
	
	int N = A.rows();
	Eigen::VectorXd dv = A.diagonal().array().square();

	Eigen::VectorXi d(1); d << 0; // represents the main diagonal

	Eigen::SparseMatrix<double> Dm = spdiags<double>(dv.cwiseInverse(),d,N,N); // \Blue{$\VD^{-\nicefrac{1}{2}}$}
	Eigen::SparseMatrix<double> Dp = spdiags<double>(dv,d,N,N);     // \Blue{$\VD^{\nicefrac{1}{2}}$%

	// Build (densely populated !) matrix \Blue{$\wh{\VA}$}
	Eigen::VectorXd c = Dp*Eigen::VectorXd::Constant(N,1); 
	Eigen::SparseMatrix<double> Ahtemp = Dm*A*Dm;
	Eigen::MatrixXd Ah = Eigen::MatrixXd(Ahtemp) + 2*c*c.transpose();


	// Compute eigenvalues; grossly inefficient \Red{!}
	Eigen::EigenSolver<Eigen::MatrixXd> solver(Ah, true);
	auto E = solver.eigenvectors().real();
	auto W = solver.eigenvalues().real();


	// Obtain eigenvector \Blue{$\Vx^{\ast}$} belonging to 2nd smallest generalized
	int maxIndex;
	W.maxCoeff(&maxIndex);

	// eigenvalue of \Blue{$\VA$} and \Blue{$\VD$}
	Eigen::VectorXd x = W.col(maxIndex);
	x = Dm*x;

	// Extract segmented image
	//Eigen::MatrixXd xs = reshape(x,m,n); 
}
