#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "copy.hpp"

#include <Eigen/Dense>

using namespace Eigen;

struct TestData {
  TestData() {
    Mesh.resize(11);
    Mesh << 0.0, 1e-6, 1e-4, 0.01, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0;

    sqrt = [](double x) { return std::sqrt(x); };

    sin = [](double x) { return std::sin(10 * x); };
  }

  VectorXd Mesh;
  std::function<double(double)> sqrt;
  std::function<double(double)> sin;
};

TestData data;

TEST_SUITE("GradedMeshes") {
  TEST_CASE("double pwlintpMaxError" *
            doctest::description("pwlintpMaxError()")) {
    double sol_sqrt = pwlintpMaxError(data.sqrt, data.Mesh);
    double stud_sqrt = pwlintpMaxError_TEST(data.sqrt, data.Mesh);

    double sol_sin = pwlintpMaxError(data.sin, data.Mesh);
    double stud_sin = pwlintpMaxError_TEST(data.sin, data.Mesh);

    CHECK(std::abs(sol_sqrt - stud_sqrt) == doctest::Approx(0.).epsilon(1e-6));
    CHECK(std::abs(sol_sin - stud_sin) == doctest::Approx(0.).epsilon(1e-6));
  }

  TEST_CASE("VectorXd polyfit" * doctest::description("polyfit helper")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("void cvgplotEquidistantMesh" *
            doctest::description("Equidistant mesh convergence plot")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("Eigen::VectorXd cvgrateEquidistantMesh" *
            doctest::description("Equidistant mesh convergence")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("void testcvgEquidistantMesh" *
            doctest::description("Equidistant mesh convergence tabulation")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("Eigen::MatrixXd cvgrateGradedMesh" *
            doctest::description("Graded mesh convergence")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }

  TEST_CASE("void testcvgGradedMesh" *
            doctest::description("Graded mesh convergence tabulation")) {
    MESSAGE("This function wasn't tested. Run the program to see its output.");
  }
}
