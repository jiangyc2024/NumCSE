#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <Eigen/Dense>
#include <vector>

#include "copy.hpp"
#include "doctest.h"

TEST_SUITE("CodeQuiz") {
  TEST_CASE("double myfunction_modified" *
            doctest::description("Test equivalent function")) {
    const std::vector<double> X = {0.0122, 0.3211, 0.7252, 1.3924, 1.5832};
    const std::vector<double> Y = {
        -4.40631932724292596503801178187, -1.13600267788666875468095440738,
        -0.321307800101441110030009440379, 0.331028876955146822425035679771,
        0.459448115306218873854504636256};
    for (unsigned int i = 0; i < X.size(); i++) {
      // const double sol = myfunction_modified(X[i]);
      const double stud = myfunction_modified_TEST(X[i]);

      CHECK(Y[i] == doctest::Approx(stud).epsilon(1e-15));
    }
  }

  TEST_CASE("double myfunction" * doctest::description("Unknown function") *
            doctest::skip()) {}
}
