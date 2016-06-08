MatrixXd lurec(const MatrixXd &A){
	int n = A.rows();
	MatrixXd result(n,n);
	if(n > 1){
		VectorXd fac = A.col(0).tail(n-1) / A(0,0);
		result.bottomRightCorner(n-1,n-1) = 
					lurec( A.bottomRightCorner(n-1,n-1)
							- fac * A.row(0).tail(n-1) );
	result.row(0) = A.row(0); result.col(0).tail(n-1) = fac;
	return result;
	}
	return A;
}
