//! in-situ recursive Gaussian elimination, no pivoting
//! right hand side in rightmost column of \Blue{$\VA$}: A.rightCols(1)
//! be aware: back substitution is not done in this code snippet!
void blockgs(MatrixXd &A, int n){
	if(n > 1){
		// \Red{rank-1 modification} of \Blue{$\VC$}
		A.bottomRightCorner(n-1, n) -= A.col(A.rows() - n).tail(n-1)
			* A.row(A.rows()-n).tail(n) / A(A.rows()-n,A.rows()-n);
		blockgs(A, n-1);
		A.col(A.rows()-n).tail(n-1).setZero(); // set \Vd  to 0
	}
}
