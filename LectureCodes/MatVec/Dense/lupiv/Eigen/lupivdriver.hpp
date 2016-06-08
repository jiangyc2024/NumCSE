void lupivdriver(const MatrixXd &A, MatrixXd &L, MatrixXd &U){
	MatrixXd Adec = A;
	lupiv(Adec);
	U = Adec.triangularView<Upper>();
	L.setIdentity();
	L += Adec.triangularView<StrictlyLower>();
}
