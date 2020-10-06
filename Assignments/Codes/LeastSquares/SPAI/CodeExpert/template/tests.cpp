#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

// includes for test data
#include <Eigen/Dense>

struct TestData {
	
};

TestData data;

TEST_SUITE("Sparse Approximate Inverse (SPAI)") {
	TEST_CASE("SparseMatrix<double> spai" * doctest::description("Sparse Approximate Inverse")) {
		
	}
	
	TEST_CASE("SparseMatrix<double> init_A" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("tuple_vector testSPAIPrecCG" * doctest::description("Iteration test")) {
		
	}
}

