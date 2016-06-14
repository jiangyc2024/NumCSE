//! in-situ Gaussian elimination, no pivoting
//! right hand side in rightmost column of \Blue{$\VA$}: A(:,end)
//! be aware: back substitution is not done in this code snippet!
void blockgs(MatrixXd &A){
	int n = A.rows();
	for(int i = 1; i < n; ++i){
		// \Red{rank-1 modification} of \Blue{$\VC$}
		A.bottomRightCorner(n-i,n-i+1) -= A.col(i-1).tail(n-i)
							* A.row(i-1).tail(n-i+1) / A(i-1,i-1);
		A.col(i-1).tail(n-i).setZero();	// set \Vd  to 0
	}
}
