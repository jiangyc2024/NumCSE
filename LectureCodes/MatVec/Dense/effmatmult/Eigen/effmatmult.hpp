//! This function compares the runtimes for the multiplication of a vector with a
//! rank-1 matrix $\Va\Vb^{\top}$, $\Va,\Vb\in\bbR^n$ using different associative
//! evaluations measurements consider minimal time for several (\texttt{nruns}) runs
MatrixXd dottenstiming(){
  int nruns = 3, minExp = 2, maxExp = 13;
  MatrixXd timings(maxExp-minExp+1,3);	// Matrix for storing recorded runtimes
  for(int i = 0; i <= maxExp-minExp; ++i){
	Timer tfool, tsmart;	// timer class
	int n = std::pow(2, minExp + i);
    VectorXd a = VectorXd::LinSpaced(n,1,n), b = VectorXd::LinSpaced(n,1,n).reverse(), 
    x = VectorXd::Random(n,1), y(n);
    for(int j = 0; j < nruns; ++j){
		tfool.start();	y = (a*b.transpose()) * x;	tfool.stop();// foolish implementation
		tsmart.start();	y = a * b.dot(x);	tsmart.stop();		// smart implementation
	}
	timings(i,0)=n; timings(i,1)=tsmart.min(); timings(i,2)=tfool.min();
  }
  return timings;
}
