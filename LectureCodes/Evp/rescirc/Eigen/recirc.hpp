#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <iostream>
#include <figure/figure.hpp>

// Ex.~\ref{ex:networlewp}: Numerical nodal analysis of the resonant circuit
// $R$, $L$, $C$ $\hat{=}$ network component parameters
void rescirc(double R, double L, double C)
{
	double Z = 1/R, K = 1/L;
	std::complex<double> i = std::complex<double>(0,1);

	// Matrices \Blue{$\VW$}, \Blue{$\VC$}, \Blue{$\VS$} for nodal analysis of circuit
	// we use matrices with std::complex<double> entries
	Eigen::Matrix3cd Wmat(3,3), Cmat(3,3), Smat(3,3);
   	Wmat <<  0, 0, 0,
		 	 0, Z, 0,
			 0, 0, 0;

   	Cmat <<	C, 0, 0,
	   		0, C, 0,
		   	0, 0, C;

   	Smat <<	 K,  -K,  0,
			-K, 2*K, -K,
		   	 0,	 -K,  K;

	// System matrix from nodal analysis
	auto Amat = [&](double w) {
		Eigen::Matrix3cd A(3,3); 
		A = Wmat+i*w*Cmat+Smat/(i*w);
		return A;
	};

	// rhs vector
	Eigen::Vector3cd b(3);
	b << C, 0, 0;

	// Scanning source currents
	int n = 199;
	// store the 3 nodal potentials for each omega (w)
	std::vector<Eigen::VectorXd> resU(3, Eigen::VectorXd(n)); 
	Eigen::VectorXd w = Eigen::VectorXd::LinSpaced(n,0.01,2);

	for (int i=0; i<n; ++i)
	{
		Eigen::Vector3cd u = Amat(w(i)).lu().solve(b);
		for (int j=0; j<3;++j) resU[j](i) = std::abs(u(j));
	}

	// plot potentials resU[i], i=0,1,2 
	mgl::Figure fig;

	for (int i=0; i<3; ++i)
	{
		std::string label = "U" + std::to_string(i+1);
		fig.plot(w, resU[i]).label(label);
	}

	fig.legend(0, 0.1);
	fig.title("resonant circuit");
	fig.xlabel("angular frequency omega of source voltage U");
	fig.ylabel("maximum nodal potential");
  	fig.setFontSize(5);
	fig.grid(false);
	fig.save("rescircpot");
	
	
	// Solving generalized eigenvalue problem \eqref{eq:nwevp1}
	Eigen::Matrix3cd Zmat = Eigen::Matrix3cd::Zero();
	Eigen::Matrix3cd Imat = Eigen::Matrix3cd::Identity();

	// Assemble \Blue{$6\times 6$}-matrices \Blue{$\VM$} and \Blue{$\VB$}
	Eigen::MatrixXcd Mmat(6,6), Bmat(6,6);
   	Mmat << Wmat, Smat, 
		 	Imat, Zmat;

	Bmat << -i*Cmat,   Zmat,
			   Zmat, i*Imat;

	// Solve \emph{generalized eigenvalue problem}, \emph{cf.} \eqref{eq:nevpgen}
	
	// Remark: Eigen doesnt support generalized complex eigenvalue problems.
	// Therefore we convert it to a standard eigenvalue problem.
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver;
	solver.compute(Bmat.inverse()*Mmat);
	auto omega = solver.eigenvalues();

	// plot frequencies
	mgl::Figure freq;
	freq.plot(omega.real(), omega.imag()," r*").label("w");
	freq.legend(0, 0.1);
	freq.title("resonances");
	freq.xlabel("Re(w)");
	freq.ylabel("Im(w)");
  	freq.setFontSize(5);
	freq.save("rescircomega");
}
