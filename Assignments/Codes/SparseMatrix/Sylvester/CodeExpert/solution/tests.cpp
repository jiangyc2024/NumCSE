#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>
#include <Eigen/Sparse>

Eigen::MatrixXd sparseSPD(int n) {
	double rate = 1.0 / n;
	// R is a random sparse matrix.
	ArrayXXd R = ArrayXXd::Random(n,n);
	R = R * (R.abs() < rate).cast<double>();
	// B is diagonally dominant.
	MatrixXd B = R.matrix();
	B.diagonal() = R.abs().matrix().rowwise().sum();
	// A is symmetric and strictly diagonally dominant.
	MatrixXd A = B + B.transpose() + MatrixXd::Identity(n,n);
	return A;
}

struct TestData {
	TestData() {
		n = 20;
		A = sparseSPD(n);
		As = A.sparseView();
		As.makeCompressed();
	}
	
	int n;
	Eigen::MatrixXd A;
	Eigen::SparseMatrix<double> As;
};

TestData data;

TEST_SUITE("Sylvester") {
	
	TEST_CASE("SparseMatrix<double> solveDiagSylvesterEq" * doctest::description("solveDiagSylvesterEq()")) {
		Eigen::ArrayXd diagA = Eigen::ArrayXd::LinSpaced(data.n, 1, data.n);
		
		Eigen::SparseMatrix<double> Xdiag_sol = solveDiagSylvesterEq(diagA);
		Eigen::SparseMatrix<double> Xdiag_stud = solveDiagSylvesterEq_TEST(diagA);
		
		CHECK((Xdiag_sol - Xdiag_stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("SparseMatrix<double> sparseKron" * doctest::description("sparseKron()")) {
		Eigen::SparseMatrix<double> sol = sparseKron(data.As);
		Eigen::SparseMatrix<double> stud = sparseKron_TEST(data.As);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("MatrixXd solveSpecialSylvesterEq" * doctest::description("solveSpecialSylvesterEq")) {
		Eigen::MatrixXd sol = solveSpecialSylvesterEq(data.As);
		Eigen::MatrixXd stud = solveSpecialSylvesterEq_TEST(data.As);
		
		CHECK((sol - stud).norm() == doctest::Approx(0.).epsilon(1e-6));
	}
	
}

