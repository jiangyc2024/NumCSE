
# include <iostream>
# include <Eigen/Dense>
# include <mgl2/mgl.h>




void TriPlot(Eigen::Matrix<double, -1, -1, Eigen::RowMajor> &T, Eigen::VectorXd &x, Eigen::VectorXd &y){
	mglGraph gr;
	mglData xd(x.data(), x.size()),
			yd(y.data(), y.size()),
			Td(T.rows(), T.cols(), T.data());
	gr.SetRanges(0,1,0,1);
	//gr.Grid("","h="); 
	gr.TriPlot(Td, xd, yd,"#b");
	gr.Plot(xd, yd, " r*");
	gr.Label(xd, yd,"%n"); // give each point an individual number
	gr.Axis();
	gr.WriteEPS("meshplot_cpp.eps");
	
}
void TriPlot(Eigen::MatrixXd &T, Eigen::VectorXd &x, Eigen::VectorXd &y){
	// NOTE: MathGL's mglData constructor needs the matrix in RowMajor!
	Eigen::Matrix<double, -1, -1, Eigen::RowMajor> TRow(T);
	TriPlot(TRow,x,y);
}


int main()
{
  	Eigen::VectorXd x(10), y(10);
	x << 1.0,0.60,0.12,0.81,0.63,0.09,0.27,0.54,0.95,0.96;
	y << 0.15,0.97,0.95,0.48,0.80,0.14,0.42,0.91,0.79,0.95;
	
	// specify triangles through indices of their vertices
	Eigen::MatrixXd T(11,3);
	T << 7, 1, 2,   5, 6, 2,    4, 1, 7,    6, 7, 2,
		6, 4, 7,   6, 5, 0,    3, 6, 0,    8, 4, 3, 
		3, 4, 6,   8, 1, 4,    9, 1, 8;

	TriPlot(T, x, y);
	
	return 0;
}
