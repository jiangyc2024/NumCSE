#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

struct TestData {
	TestData() {
		f = [] (double x) {
			return std::log(std::abs(x) + 1) * std::pow(x, 3);
		};
		
		g = [] (double x, double y) {
			return x * (y - 1) * std::sin(6 * x * y);
		};
		
		Q.nodes.resize(5);
		Q.nodes << -1., -std::sqrt(3./7.), 0, std::sqrt(3./7.), 1.;
		Q.weights.resize(5);
		Q.weights << 0.1, 49./90., 32./45., 49./90., 0.1;
		
		a = -2;
		b = 3;
		N = 5;
	}
	
	std::function<double (double)> f;
	std::function<double (double, double)> g;
	
    QuadRule Q;
    
    double a;
    double b;
    int N;
};

TestData data;

TEST_SUITE("NestedQuad") {
	TEST_CASE("double evalquad" * doctest::description("evalquad()")) {
		double sol = evalquad(data.a, data.b, data.f, data.Q);
		double stud = evalquad_TEST(data.a, data.b, data.f, data.Q);
		
		CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("double gaussquadtriangle" * doctest::description("evalquad()")) {
		double sol = gaussquadtriangle(data.g, data.N);
		double stud = gaussquadtriangle_TEST(data.g, data.N);
		
		CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
	}
	
	TEST_CASE("void convtest2DQuad" * doctest::description("convtest2DQuad() for 2 nodes") * doctest::skip()) {}
}

