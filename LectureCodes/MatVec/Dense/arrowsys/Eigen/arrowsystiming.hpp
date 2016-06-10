MatrixXd arrowsystiming(){
	std::vector<int> n = {8,16,32,64,128,256,512,1024,2048,4096};
	int nruns = 3;
	MatrixXd times(n.size(),6);	
	for(int i = 0; i < n.size(); ++i){
		Timer t1, t2, t3, t4;	// timer class
		double alpha = 2;
		VectorXd b = VectorXd::Ones(n[i],1);
		VectorXd c = VectorXd::LinSpaced(n[i],1,n[i]);
		VectorXd d = -b;
		VectorXd y = VectorXd::Constant(n[i]+1,-1).binaryExpr(
				VectorXd::LinSpaced(n[i]+1,1,n[i]+1), 
				[](double x, double y){return pow(x,y);} ).array();
		VectorXd x1(n[i]+1), x2(n[i]+1), x3(n[i]+1), x4(n[i]+1);
		for(int j = 0; j < nruns; ++j){
			t1.start();	x1 = arrowsys_slow(d,c,b,alpha,y);	t2.stop();
			t2.start();	x2 = arrowsys_fast(d,c,b,alpha,y);	t2.stop();
			t3.start();
			x3 = arrowsys_sparse<SparseLU<SparseMatrix<double> > >(d,c,b,alpha,y);
			t3.stop();
			t4.start();
			x4 = arrowsys_sparse<BiCGSTAB<SparseMatrix<double> > >(d,c,b,alpha,y);
			t4.stop();
		}
		times(i,0)=n[i]; times(i,1)=t1.min(); times(i,2)=t2.min();
		times(i,3)=t3.min();times(i,4)=t4.min();times(i,5)=(x4-x3).norm();
	}
	return times;
}
