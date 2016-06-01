//! Eigen script for timing numerical solution of linear systems
MatrixXd gausstiming(){
	std::vector<int> n = {8,16,32,64,128,256,512,1024,2048,4096,8192};
	int nruns = 3;
	MatrixXd times(n.size(),3);
	Timer timer;	// timer class
	for(int i = 0; i < n.size(); ++i){
		MatrixXd A = MatrixXd::Random(n[i],n[i]) + n[i]*MatrixXd::Identity(n[i],n[i]);
		VectorXd b = VectorXd::Random(n[i]);
		VectorXd x(n[i]);
		double t1 = std::numeric_limits<double>::max(); double t2 = t1;
		for(int j = 0; j < nruns; ++j){
			timer.start();
			x = A.lu().solve(b);	// Eigen implementation
			timer.stop();
			t1 = std::min(t1, timer.duration());
			#ifndef EIGEN_USE_MKL_ALL
			if(n[i] <= 4096){		// Prevent long runs
			timer.start();
			gausselimsolve(A,b,x);	// own gauss elimination
			timer.stop();
			t2 = std::min(t2, timer.duration());
			}else t2=0;				// is omitted in plot
			#else
			t2 = 0;
			#endif
		}
		times(i,0) = n[i]; times(i,1) = t1; times(i,2) = t2;
	}
	return times;
}
