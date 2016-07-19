#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/KroneckerProduct> 


/* Computes vibration modes of a truss structure, see Ex.~\ref{ex:truss}. Mass point
 positions passed in the $n\times 2$-matrix \texttt{poss} and the connectivity encoded in
 the sparse symmetric matrix \texttt{top}. In addition \texttt{top(i,j)} also stores the
 Young's moduli \Blue{$\alpha_{ij}$}.
 The \Blue{$2n$} resonant frequencies are returned in the vector \texttt{lambda}, the
 eigenmodes in the column of \texttt{V}, where entries at odd positions contain the
 \Blue{$x_{1}$}-coordinates, entries at even positions the \Blue{$x_{2}$}-coordinates
*/
Eigen::EigenSolver<Eigen::MatrixXd> trussvib(const Eigen::MatrixXd& pos, const Eigen::MatrixXd& con)
{

	int n = pos.rows(); // no. of point masses

	// Assembly of stiffness matrix according to \eqref{eq:stiffmat}
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*n,2*n);
	Eigen::MatrixXd top(n,n); // top(i,j) = top(j,i) = 1 iff there is a (i,j) in con

	//[Iidx,Jidx] = find(top); idx = [Iidx,Jidx]; % Find connected masses
	for (int k=0; k<con.rows(); ++k)
	{
		int i = con(k,0)-1;
		int j = con(k,1)-1;

		Eigen::VectorXd dp(2);
		dp = pos.row(j) - pos.row(i); // \Blue{$\Delta\Vp^{ji}$}
		double lij = dp.norm();       // \Blue{$l_{ij}$} 

		A.block(2*i, 2*j, 2, 2) = -(dp*dp.transpose())/std::pow(lij,3);
		A.block(2*j, 2*i, 2, 2) = -(dp*dp.transpose())/std::pow(lij,3);

		top(i,j) = 1;
		top(j,i) = 1;
	}

	// Set Young's moduli \Blue{$\alpha_{ij}$} (stored in \texttt{top} matrix)
	Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(2,2);
	A = A.cwiseProduct(Eigen::kroneckerProduct(top, ones).eval());

	//Set $2\times2$ diagonal blocks
	for (int i=0; i<2*n; ++i)
	{
		double sumL = 0; // left value of block matrix
		double sumR = 0; // right value of block matrix

		for (int j=0; j<2*n; ++j)
		{
			if (j%2==0) sumL += A(i, j);
			else sumR += A(i,j);
		}

		if (i%2==0) // write upper part of 2x2 block
		{
			A(i,i) = -sumL;
			A(i,i+1) = -sumR;
		}
		else // write lower part of 2x2 block
		{
			A(i,i-1) = -sumL;
			A(i,i) = -sumR;
		}
	}

	// Compute eigenvalues and eigenmodes
	Eigen::EigenSolver<Eigen::MatrixXd> ev(A);
	return ev;
}


