#include <Eigen/Dense>
#include <figure.hpp>
#include <iostream>
#include <eigensolversort.hpp>

/* monitor power iteration with Ritz projection for computing
 * the two largest (in modulus) eigenvalues and associated eigenvectors
 * of a symmetric matrix with prescribed eigenvalues passed in \texttt{d} */
void sspowitrpex(const Eigen::VectorXd& d, int maxit = 20)
{
	// Generate matrix
	int n = d.rows();
	Eigen::VectorXd oneToN = Eigen::VectorXd::LinSpaced(n,1,n);
	Eigen::MatrixXd Z = Eigen::MatrixXd::Ones(n,n);
	Z += oneToN.array().sqrt().matrix().asDiagonal();

	auto qr = Z.householderQr();   // generate orthogonal matrix
	Eigen::MatrixXd Q = qr.householderQ();
	Eigen::MatrixXd A = Q*d.asDiagonal()*Q.transpose(); // ``synthetic'' \Blue{$\VA=\VA^T$} with spectrum \Blue{$\sigma(\VA) =\{d_1,\ldots,d_n\}$}

	// Compute ``exact'' eigenvectors and eigenvalues
	Eigen::EigenSolver<Eigen::MatrixXd> solver(A);
	Eigen::VectorXd D; // eigenvalues
	Eigen::MatrixXd V; // eigenvectors
	std::tie(D, V) = eigensolversort(solver);

	Eigen::VectorXd v_ex = V.col(n-1); // eigenvector belonging to largest eigenvalue
	Eigen::VectorXd w_ex = V.col(n-2); // eigenvector belonging to second largest eigenvalue
	double lv_ex = D(n-1); // largest eigenvalue
	double lw_ex = D(n-2); // second largest eigenvalue

	// (Arbitrary) initial guess for eigenvectors
	Eigen::VectorXd v = Eigen::VectorXd::Ones(n);
	Eigen::VectorXd w = -1*Eigen::VectorXd::Ones(n);
	v.normalize(); w.normalize();

	Eigen::MatrixXd result(maxit, 5);
	Eigen::VectorXd v_new;
	Eigen::VectorXd w_new;


	for (int k=0; k <maxit; ++k)
	{
		// ``power iteration'', \emph{cf.} \eqref{eq:dirpotmeth}
		v_new = A*v; w_new = A*w; 

		// orthogonalization, \emph{cf.} Rem.~\ref{rem:QRorth}
		Eigen::MatrixXd B = Eigen::MatrixXd::Zero(v.rows(), 2);
		B.col(0) = v_new;
		B.col(1) = w_new;

		// Solve Ritz projected eigenvalue problem
		auto qr = B.householderQr();
		Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n,2); // economy size qr

		Eigen::EigenSolver<Eigen::MatrixXd> esolver(Q.transpose()*A*Q);

		// recover approximate eigenvectors
		Eigen::VectorXd D; // eigenvalues
		Eigen::MatrixXd U; // eigenvectors
		std::tie(D, U) = eigensolversort(esolver, true);
		v = Q*U.col(1);
		w = Q*U.col(0);
		
		// Record errors in eigenvalue and eigenvector approximations. Note that the 
		// direction of the eigenvectors is not specified.
		result.row(k) << 
				k, 
				std::abs(D(1) - lv_ex),
				std::abs(D(0) - lw_ex),
				std::min((v-v_ex).norm(), (v+v_ex).norm()), 
				std::min((w-w_ex).norm(), (w+w_ex).norm());
	}

	std::cout << result << std::endl;

	// plot errors
	mgl::Figure fig1;
  	fig1.setlog(false, true);
	fig1.plot(result.col(0), result.col(1), " m-+").label("error in lambda_n");
	fig1.plot(result.col(0), result.col(2), " r-+").label("error in lambda_n-1");
	fig1.plot(result.col(0), result.col(3), " k-*").label("error in v");
	fig1.plot(result.col(0), result.col(4), " b-*").label("error in w");
	fig1.title("sspowitrp");
	fig1.legend(1,1);
	fig1.xlabel("power iteration step");
	fig1.ylabel("error");
	fig1.setFontSize(5);
	fig1.save("sspowitrpcvg");

	// plot rates
	Eigen::MatrixXd rates = result.bottomRightCorner(maxit-1,4).cwiseQuotient(result.topRightCorner(maxit-1,4));
	mgl::Figure fig2;
	fig2.plot(result.col(0).tail(maxit-1), rates.col(0), " m-+").label("error in lambda_n");
	fig2.plot(result.col(0).tail(maxit-1), rates.col(1), " r-+").label("error in lambda_n-1");
	fig2.plot(result.col(0).tail(maxit-1), rates.col(2), " k-*").label("error v");
	fig2.plot(result.col(0).tail(maxit-1), rates.col(3), " b-*").label("error w");
	fig2.xlabel("power iteration step");
	fig2.ylabel("error quotient");
	fig2.title("sspowit rates");
	fig2.legend(1,1);
	fig2.setFontSize(5);
	fig2.save("sspowitrpcvgrates");
} 
