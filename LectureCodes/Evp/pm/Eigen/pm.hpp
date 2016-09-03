#include <figure/figure.hpp>
#include <Eigen/Dense>


// Demonstration of direct power method for Ex.~\ref{ex:pm}
void directPowerMethod()
{
	int maxit = 30;
	int n = 10;
	Eigen::VectorXd d = Eigen::VectorXd::LinSpaced(n,1,n);
	
	// Initialize the matrix \Blue{$\VA$}
	
	Eigen::MatrixXd S = Eigen::MatrixXd::Ones(n,n);
   	S += d.reverse().asDiagonal();
	S = S.triangularView<Eigen::Upper>();
	Eigen::MatrixXd A = S*d.asDiagonal()*S.inverse();


	// This calculates the exact eigenvalues (for error calculation)
	Eigen::EigenSolver<Eigen::MatrixXd> solver(A);
	Eigen::VectorXd ev = solver.eigenvectors().col(n-1).real(); ev.normalize(); // normalized eigenvector of max eigenvalue
	double ew = d(n-1); // max eigenvalue
	int sgn = 1; // because ew > 0


	Eigen::VectorXd z = Eigen::VectorXd::Random(n,1);
	z.normalize();
	int s = 1;

	Eigen::MatrixXd res(maxit,5);

	// Actual direct power iteration
	for (int i=0; i < maxit; ++i)
	{
		Eigen::VectorXd w = A*z;
		double l = w.norm();
		double rq = w.dot(z);
		z = w/l;
		res.row(i) << i,l,(s*z-ev).norm(), std::abs(l-std::abs(ew)), std::abs(sgn*rq-ew);
		s = s*sgn;
	}


	// Plot the result
	
	mgl::Figure fig;
	fig.plot(res.col(0), res.col(2)).label("(s*z-ev).norm()");
	fig.plot(res.col(0), res.col(3)).label("abs(l-abs(ew))");
	fig.plot(res.col(0), res.col(4)).label("abs(sign*(dot(w,z)-ew)");

	fig.legend(1, 1);
  	fig.setlog(false, true, false);
	fig.title("direct power method");
	fig.xlabel("iteration step k");
	fig.ylabel("errors");
	fig.setFontSize(5);
	fig.save("pm1");
}
