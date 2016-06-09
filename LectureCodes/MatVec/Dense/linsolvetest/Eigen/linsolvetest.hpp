//! Eigen script: assessing the gain from using special properties
//! of system matrices in Eigen
MatrixXd gausstiming(){
	std::vector<int> n = {16,32,64,128,256,512,1024,2048,4096,8192};
	int nruns = 3;
	MatrixXd times(n.size(),3);
	Timer timer;	// timer class
	for(int i = 0; i < n.size(); ++i){
		MatrixXd A = VectorXd::LinSpaced(n[i],1,n[i]).asDiagonal() ;
		A += MatrixXd::Ones(n[i],n[i]).triangularView<Upper>();
		VectorXd b = VectorXd::Random(n[i]);
		VectorXd x1(n[i]), x2(n[i]);
		double t1 = std::numeric_limits<double>::max(); double t2 = t1;
		for(int j = 0; j < nruns; ++j){
			timer.start();
			x1 = A.lu().solve(b);
			timer.stop();
			t1 = std::min(t1, timer.duration());
			timer.start();
			x2 = A.triangularView<Upper>().solve(b);
			timer.stop();
			t2 = std::min(t2, timer.duration());
		}
		times(i,0) = n[i]; times(i,1) = t1; times(i,2) = t2;
	}
	return times;
}
