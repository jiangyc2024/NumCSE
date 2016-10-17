#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/LU>

using namespace Eigen;

void spai(const MatrixXd& A, MatrixXd& B) { // $A$ must be square...
	B = A;
	for(size_t i=0; i<A.cols(); ++i) {
		std::vector<size_t> idx;
		MatrixXd tempMat;
		for(size_t j=0; j<A.rows(); ++j) {
			if(A(j,i) != 0) {
				idx.push_back(j);
				tempMat.conservativeResize(A.rows(), tempMat.cols()+1);
				tempMat.col(tempMat.cols()-1) = A.col(j);
			}
		}
		VectorXd tempVec = (tempMat.transpose() * tempMat).fullPivLu().solve(tempMat.row(i));
		for(size_t j=0; j<idx.size(); ++j) {
			B(j,i) = tempVec(j);
		}
	}
}

int main() {

}
