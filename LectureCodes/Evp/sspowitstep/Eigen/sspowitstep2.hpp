#include <Eigen/Dense>

inline void sspowitstep2(Eigen::VectorXd &v, Eigen::VectorXd &w)
{
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(v.rows(), 2);
	A.col(0) = v;
	A.col(1) = w;
	auto qr = A.householderQr();
	
	Eigen::MatrixXd Q = qr.householderQ();
	
	v = Q.col(0);
	w = Q.col(1);
}
