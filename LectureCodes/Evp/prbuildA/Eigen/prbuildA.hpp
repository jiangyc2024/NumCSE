#include <Eigen/Dense>


Eigen::MatrixXd prbuildA(const Eigen::MatrixXi &G, double d) 
{
	int N = G.rows();

	Eigen::VectorXi sum = G.colwise().sum();

	Eigen::VectorXd s = Eigen::VectorXd::Zero(N);
	Eigen::VectorXd ds = Eigen::VectorXd::Ones(N) / N;

	for (int i=0; i<N; ++i)
	{
		if (sum(i) > 0)
		{
			ds(i) = d/N;
			s(i) = d/N;
		}
	}

	Eigen::MatrixXd A = Eigen::MatrixXd::Ones(N,N);
	A *= ds.asDiagonal();
	A += (1-d)*G.cast<double>()*s.asDiagonal();

	return A;
}
