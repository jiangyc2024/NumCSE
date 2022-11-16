#include <Eigen/Dense>

/* Computes one step of subspace power iteration. A needs to be symmetric. 
 * The arguments v and w will be modified.
*/
inline void sspowitstep(const Eigen::MatrixXd &A, Eigen::VectorXd &v, Eigen::VectorXd &w)
{
	v = A*v;
	w = A*w;
	v.normalize();
	w = w - v.dot(w)*v;
	w.normalize();
}
