#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

int main () {
	MatrixXd A(3,3);
	A.setRandom();
	/* SAM_LISTING_BEGIN_0 */
	Eigen::PartialPivLU<MatrixXd> lu(A);
	MatrixXd L = MatrixXd::Identity(3,3);
	L.triangularView<StrictlyLower>() += lu.matrixLU();
	MatrixXd U = lu.matrixLU().triangularView<Upper>();
	MatrixXd P = lu.permutationP();
	/* SAM_LISTING_END_0 */
	std::cout << A << std::endl;
	std::cout << L << std::endl;
	std::cout << U << std::endl;
	std::cout << P*L*U-A << std::endl;
	
	return 0;
}
