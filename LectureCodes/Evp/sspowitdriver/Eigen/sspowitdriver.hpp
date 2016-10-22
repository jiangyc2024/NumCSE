#include <Eigen/Dense>
#include <Eigen/QR>
#include <figure.hpp>
#include <eigensolversort.hpp>


/* monitor power iteration with orthogonal projection for finding
 * the two largest (in modulus) eigenvalues and associated eigenvectors
 * of a symmetric matrix with prescribed eigenvalues passed in \texttt{d} */
void ssppowitdriver(const Eigen::VectorXd &d, int maxit = 20)
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
	double lv, lw;
	for (int k=0; k <maxit; ++k)
	{
		v_new = A*v; w_new = A*w; // ``power iteration'', \emph{cf.} \eqref{eq:dirpotmeth}

		// Rayleigh quotients provide approximate eigenvalues
		lv = v_new.dot(v); lw = w_new.dot(w); 

		// orthogonalization, \emph{cf.} Gram-Schmidt orthogonalization \eqref{GRS}: \Blue{$\Vw\perp \Vv$}
		v = v_new/v_new.norm(); w = w_new - v.dot(w_new)*v; w.normalize();

		// Record errors in eigenvalue and eigenvector approximations. Note that the 
		// direction of the eigenvectors is not specified.
		result.row(k) << k, 
				std::abs(lv-lv_ex),
				std::abs(lw-lw_ex), 
				std::min((v-v_ex).norm(), (v+v_ex).norm()), 
				std::min((w-w_ex).norm(), (w+w_ex).norm());
	}

	// plot errors
	mgl::Figure fig1;
  	fig1.setlog(false, true);
	fig1.plot(result.col(0), result.col(1), " m-+").label("error in lambda_n");
	fig1.plot(result.col(0), result.col(2), " r-+").label("error in lambda_n-1");
	fig1.plot(result.col(0), result.col(3), " k-*").label("error in v");
	fig1.plot(result.col(0), result.col(4), " b-*").label("error in w");
	fig1.title("sspowit");
	fig1.legend(1,1);
	fig1.xlabel("power iteration step");
	fig1.ylabel("error");
	fig1.setFontSize(5);
	fig1.save("sspowitcvg1");

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
	fig2.save("sspowitcvgrates1");
}


     
