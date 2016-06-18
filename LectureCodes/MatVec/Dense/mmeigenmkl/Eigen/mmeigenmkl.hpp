//! script for timing different implementations of matrix multiplications
void mmeigenmkl(){
  int nruns = 3, minExp = 6, maxExp = 13;
  MatrixXd timings(maxExp-minExp+1,2);
  for(int p = 0; p <= maxExp-minExp; ++p){
	Timer t1;	// timer class
	int n = std::pow(2, minExp + p);
    MatrixXd A = MatrixXd::Random(n,n);
    MatrixXd B = MatrixXd::Random(n,n);
    MatrixXd C = MatrixXd::Zero(n,n);
    for(int q = 0; q < nruns; ++q){
		t1.start();
		C = A * B;
		t1.stop();
	}
	timings(p,0)=n; timings(p,1)=t1.min();
  }
  std::cout << timings << std::endl;
}
