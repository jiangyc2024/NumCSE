MatrixXd arrowsystiming(){
	std::vector<int> n = {8,16,32,64,128,256,512,1024,2048,4096};
	int nruns = 3;
	MatrixXd times(n.size(),6);
	Timer timer;	// timer class
	for(int i = 0; i < n.size(); ++i){
		double alpha = 2;
		VectorXd b = VectorXd::Ones(n[i],1);
		VectorXd c = VectorXd::LinSpaced(n[i],1,n[i]);
		VectorXd d = -b;
		VectorXd y = VectorXd::Constant(n[i]+1,-1).binaryExpr(
				VectorXd::LinSpaced(n[i]+1,1,n[i]+1), 
				[](double x, double y){return pow(x,y);} ).array();
		VectorXd x1(n[i]+1), x2(n[i]+1), x3(n[i]+1), x4(n[i]+1);
		double t1=numeric_limits<double>::max(), t2=t1, t3=t1, t4=t1;
		for(int j = 0; j < nruns; ++j){
			timer.start();
			x1 = arrowsys_slow(d,c,b,alpha,y);
			timer.stop();
			t1 = std::min(t1, timer.duration());
			timer.start();
			x2 = arrowsys_fast(d,c,b,alpha,y);
			timer.stop();
			t2 = std::min(t2, timer.duration());
			timer.start();
			x3 = arrowsys_sparse<SparseLU<SparseMatrix<double> > >(d,c,b,alpha,y);
			timer.stop();
			t3 = std::min(t3, timer.duration());
			timer.start();
			x4 = arrowsys_sparse<BiCGSTAB<SparseMatrix<double> > >(d,c,b,alpha,y);
			timer.stop();
			t4 = std::min(t4, timer.duration());
		}
		times(i,0)=n[i]; times(i,1)=t1; times(i,2)=t2; times(i,3)=t3;
		times(i,4)=t4;	times(i,5)=(x4-x3).norm(); // error norm
	}
	return times;
}
