//! in-situ recursive Gaussian elimination, no pivoting
//! right hand side in rightmost column of \Blue{$\VA$}: A.rightCols(1)
//! be aware: back substitution is not done in this code snippet!
void blockgs(MatrixXd &A, int n){
	if(n > 1){
		// \Red{rank-1 modification} of \Blue{$\VC$}
		A.bottomRightCorner(n-1, n-1) -= A.block(A.rows()-n+1, A.cols()-n, n-1, 1)	* (A.block(A.rows()-n,A.cols()-n+1,1,n-1)/A(A.rows()-n,A.cols()-n));
		blockgs(A, n-1);
		// set \Vd  to 0
		A.block(A.rows()-n+1,A.rows()-n,n-1,1) = MatrixXd::Zero(n-1,1);
	}
}
