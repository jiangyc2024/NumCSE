//! in-situ Gaussian elimination, no pivoting
//! be aware: back substitution is not done in this code snippet!
void blockgs(MatrixXd &A){
	int n = A.rows();
	for(int i = 1; i < n; ++i){
		// \Red{rank-1 modification} of \Blue{$\VC$}
		A.bottomRightCorner(n-i,n-i) -= A.block(i,i-1,n-i,1) 
							* (A.block(i-1,i,1,n-i)/A(i-1,i-1));
		A.block(i,i-1, n-i,1) = MatrixXd::Zero(n-i,1);	// set \Vd  to 0
	}
}
