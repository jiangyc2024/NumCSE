//! post-processing: extract \Blue{$\VL$} and \Blue{$\VU$}
void lurecdriver(const MatrixXd &A, MatrixXd &L, MatrixXd &U){
	MatrixXd A_dec = lurec(A);
	// post-processing: 
	//extract \Blue{$\VL$} and \Blue{$\VU$}
	U = A_dec.triangularView<Upper>();
	L.setIdentity();
	L += A_dec.triangularView<StrictlyLower>();
}
