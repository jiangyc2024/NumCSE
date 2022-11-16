#include "imgsegmat.hpp"
#include "imread.hpp"
#include "spdiags.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

namespace seg {


// segmentation of P: build matrices, see \cref{mc:imgsegmat} and \eqref{eq:segmat}
inline void seg(Eigen::MatrixXd &P)
{
	Eigen::SparseMatrix<double> D;
	Eigen::SparseMatrix<double> A;

	auto Sfun = [](double x, double y) 
		{ return std::exp(-0.1*std::pow(x-y, 2)); };

	std::tie(A,D) = imgsegmat::imgsegmat(P, Sfun);

	// Build scaling matrics 
	
	const Eigen::Index N = A.rows();
	const Eigen::VectorXd dv = A.diagonal().array().square();

	Eigen::VectorXi d(1); d << 0; // represents the main diagonal

	const Eigen::SparseMatrix<double> Dm = spdiags<double>(dv.cwiseInverse(),d,N,N); // \Blue{$\VD^{-\nicefrac{1}{2}}$}
	const Eigen::SparseMatrix<double> Dp = spdiags<double>(dv,d,N,N);     // \Blue{$\VD^{\nicefrac{1}{2}}$%

	// Build (densely populated !) matrix \Blue{$\wh{\VA}$}
	const Eigen::VectorXd c = Dp*Eigen::VectorXd::Constant(N,1); 
	const Eigen::SparseMatrix<double> Ahtemp = Dm*A*Dm;
	const Eigen::MatrixXd Ah = Eigen::MatrixXd(Ahtemp) + 2*c*c.transpose();


	// Compute eigenvalues; grossly inefficient \Red{!}
	const Eigen::EigenSolver<Eigen::MatrixXd> solver(Ah, true);
	auto W = solver.eigenvalues().real();


	// Obtain eigenvector \Blue{$\Vx^{\ast}$} belonging to 2nd smallest generalized
	int maxIndex = -1;
	W.maxCoeff(&maxIndex);

	// eigenvalue of \Blue{$\VA$} and \Blue{$\VD$}
	Eigen::VectorXd x = W.col(maxIndex);
	x = Dm*x;

	// Extract segmented image
	//Eigen::MatrixXd xs = reshape(x,P.rows(),P.cols()); 
}


} //namespace seg
