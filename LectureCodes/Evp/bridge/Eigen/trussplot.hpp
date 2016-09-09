#include <figure.hpp>


void trussplot(mgl::Figure& fig, const Eigen::MatrixXd& pos, const Eigen::MatrixXd& con, char style[])
{
	Eigen::VectorXd x(2);
	Eigen::VectorXd y(2);

	for (int i=0;i<con.rows(); ++i)
	{
		int k = con(i,0)-1;
		int l = con(i,1)-1;

		x << pos(k,0), pos(l, 0);
		y << pos(k,1), pos(l, 1);

		fig.plot(x,y, style);
	}

	fig.ranges(-1, 6, -1, 2);
}
