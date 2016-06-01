#pragma once
//! This function compares the runtimes for the multiplication of a vector with a
//! rank-1 matrix $\Va\Vb^{\top}$, $\Va,\Vb\in\bbR^n$ using different associative
//! evaluations measurements consider minimal time for several (\texttt{nruns}) runs
void dottenstiming(){
  int nruns = 3, minExp = 2, maxExp = 13;
  MatrixXd timings(maxExp-minExp+1,3);	// Matrix for storing recorded runtimes
  Timer timer;	// timer class
  for(int i = 0; i <= maxExp-minExp; ++i){
	int n = std::pow(2, minExp + i);
    VectorXd a = VectorXd::LinSpaced(n,1,n), b = VectorXd::LinSpaced(n,1,n).reverse(), 
    x = VectorXd::Random(n,1), y(n);
    double tfool = std::numeric_limits<double>::max();
    double tsmart = tfool;
    for(int j = 0; j < nruns; ++j){
		timer.start();
		y = (a*b.transpose()) * x;	// foolish implementation
		timer.stop();
		tfool = std::min(tfool, timer.duration());
		timer.start();
		y = a * b.dot(x);			// smart implementation
		timer.stop();
		tsmart = std::min(tsmart, timer.duration());
	}
	timings(i,0)=n; timings(i,1)=tsmart; timings(i,2)=tfool;
  }
  std::cout << timings << std::endl;
  //Plotting
  mgl::Figure fig;
  fig.setFontSize(4);
  fig.title("Timings for rank 1 matrix-vector multiplications");
  fig.setlog(true, true);
  fig.plot(timings.col(0),timings.col(1), " +b").label("efficient evaluation");
  fig.plot(timings.col(0),timings.col(2)," ^r").label("slow evaluation");
  fig.plot(timings.col(0),timings.col(0).array().pow(1).matrix()*timings(timings.rows()-1,1)
  /(std::pow(timings(timings.rows()-1,0),1)), "l;").label("O(n)");
  fig.plot(timings.col(0),timings.col(0).array().pow(2).matrix()*timings(timings.rows()-1,2)
  /(std::pow(timings(timings.rows()-1,0),2)), "h;").label("O(n^2)");
  fig.xlabel("vector length n");
  fig.ylabel("time [s]");
  fig.legend(0.05,0.95);
  fig.save("dottenstiming");
}
