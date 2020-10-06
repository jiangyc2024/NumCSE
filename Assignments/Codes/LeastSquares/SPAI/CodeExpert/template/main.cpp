#include <iostream>
#include <iomanip>

#include "spai.hpp"

int main(){
	std::cout << "Test spai():\n";
	const unsigned int n = 60;
	SparseMatrix<double> M(5, 5);
	M.coeffRef(3, 4) = 1;
	M.coeffRef(4, 4) = 2;
	M.coeffRef(1, 4) = 3;
	M.coeffRef(3, 3) = 4;
	M.coeffRef(3, 2) = 4;
	M.coeffRef(2, 3) = 4;
	M.coeffRef(2, 2) = 5;
	M.coeffRef(3, 1) = 6;
	M.coeffRef(0, 0) = 9;
	
	SparseMatrix<double> N = spai(M);
	SparseMatrix<double> I(5, 5);
	I.setIdentity();
	
	std::cout << "M = " << std::endl
	<< M
	<< std::endl;
	std::cout << "N = " << std::endl
	<< N
	<< std::endl;
	
	std::cout << "Error: "
	<< (I - M * N).norm()
	<< std::endl;
	
	SparseMatrix<double> M2(n * n, n * n);
	SparseMatrix<double> I2(n, n);
	I2.setIdentity();
	
	MatrixXd R = MatrixXd::Random(n, n);
	M2 = kroneckerProduct(R, I2);

	SparseMatrix<double> N2 = spai(M2);
	
	SparseMatrix<double> Ibig(n * n, n * n);
	Ibig.setIdentity();
	
	std::cout << "Error (n = " << n * n << "): "
	<< (Ibig - M2 * N2).norm()
	<< std::endl;
	
	tuple_vector tv = testSPAIPrecCG(5);
	std::cout << "Table for the number of iterations:\n";
	std::cout << std::setw(30) << "N" <<
	std::setw(30) << "#it.'s w/out preconditioner" <<
	std::setw(30) << "#it.'s w/ preconditioner\n";
	for (auto& tuple: tv) {
		std::cout << std::setw(30) << std::get<0>(tuple) <<
		std::setw(30) << std::get<1>(tuple) <<
		std::setw(30) << std::get<2>(tuple) <<
		std::endl;
	}
}
