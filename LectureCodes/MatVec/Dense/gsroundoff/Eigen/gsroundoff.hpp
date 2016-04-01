//! Eigen Function demonstrating the effect of roundoff on the result of Gram-Schmidt orthogonalization
//! A is 10x10 special matrix the so-called Hilbert matrix: $\MAc{i}{j} = (i+j-1)^{-1}$
void gsroundoff(MatrixXd& A){
	MatrixXd Q = gramschmidt(A); // Gram-Schmidt orthogonalization of columns of A
	// \Magenta{Test orthonormality} of column of Q, which should be an \cor{orthogonal}
	// matrix according to theory
	std::cout << std::setprecision(4) << std::fixed << "I = " 
	<< std::endl << Q*Q.transpose() << std::endl; // Should be the identity matrix, but isn't !
	// Eigens's internal Gram-Schmidt orthogonalization by \cor{QR-decomposition}
	HouseholderQR<MatrixXd> qr(A.rows(),A.cols());
	qr.compute(A); MatrixXd Q1 = qr.householderQ();
	std::cout << "I1 = " << std::endl << Q1*Q1.transpose() << std::endl; // Is identity matrix
}
