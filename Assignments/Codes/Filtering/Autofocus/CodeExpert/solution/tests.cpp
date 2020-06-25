#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

struct TestData {
	TestData() {
        // TODO: use the correct one
        M.resize(2, 2);
        M << 1, 2,
             3, 4;
	}
    Eigen::MatrixXd M;
};

TestData data;

TEST_SUITE("Autofocus") {
	TEST_CASE("MatrixXd set_focus" * doctest::description("skipped") * doctest::skip()) {}
	
	TEST_CASE("double high_frequency_content" * doctest::description("Compute V")) {		
		double sol = high_frequency_content(data.M);
		double stud = high_frequency_content_TEST(data.M);
		
        CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
	}

    TEST_CASE("double autofocus" * doctest::description("Find most focused image")) {		
		double sol = autofocus();
		double stud = autofocus_TEST();
		
		CHECK(std::abs(sol - stud) == doctest::Approx(0.).epsilon(1e-6));
	}
		
	TEST_CASE("void save_image" * doctest::description("Save image")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}

    TEST_CASE("void plot_freq" * doctest::description("Plot spectrum")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}

    TEST_CASE("void plotV" * doctest::description("Plot V")) {
		MESSAGE("This function wasn't tested. Run the program to see its output.");
	}
}

