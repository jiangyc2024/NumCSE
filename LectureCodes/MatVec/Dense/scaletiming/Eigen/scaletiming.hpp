#pragma once
//! script for timing a smart and foolish way to carry out
//! miltiplication with a scaling matrix
void scaletiming(){
  int nruns = 3, minExp = 2, maxExp = 14;
  MatrixXd timings(maxExp-minExp+1,3);
  for(int i = 0; i <= maxExp-minExp; ++i){
	Timer tbad, tgood;	// timer class
	int n = std::pow(2, minExp + i);
    VectorXd d = VectorXd::Random(n,1), x = VectorXd::Random(n,1), y(n);
    for(int j = 0; j < nruns; ++j){
		MatrixXd D = d.asDiagonal(); // d.asDiagonal()*x would be optimized
		tbad.start(); y = D*x; 	tbad.stop();	// matrix vector multiplication
		tgood.start();y=d.cwiseProduct(x);	tgood.stop(); // comp. wise mult.
	}
	timings(i,0)=n; timings(i,1)=tgood.min(); timings(i,2)=tbad.min();
  }
  std::cout << timings << std::endl;
  //Plotting
  mgl::Figure fig;
  fig.setFontSize(4);
  fig.title("Timings for different ways to do scaling");
  fig.setlog(true, true);
  fig.plot(timings.col(0),timings.col(1), " +b").label("d.*x");
  fig.plot(timings.col(0),timings.col(2)," ^r").label("asDiagonal()*x");
  fig.xlabel("vector length n");  fig.ylabel("time [s]");
  fig.legend(0.05,0.95);  fig.save("scaletiming");
}