void lurecdriver(const MatrixXd &A, 
		MatrixXd &L, MatrixXd &U){
	MatrixXd Adec = lurec(A);
	// post-processing: 
	//extract \Blue{$\VL$} and \Blue{$\VU$}
	U = Adec.triangularView<Upper>();
	L.setIdentity();
	L += Adec.triangularView<StrictlyLower>();
}
