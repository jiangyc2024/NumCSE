#include <Eigen/Dense>

void sspowitstep1(Eigen::VectorXd &v, Eigen::VectorXd &w)
{
	v.normalize();
	w = w- v.dot(w)*v;
	w.normalize();
}
