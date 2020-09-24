#include "distfitting.hpp"

int main() {
  /*
   *	run initA
   */
  unsigned int n = 4;
  MatrixXd Adense = MatrixXd(initA(n));

  std::cout << "A(4) is:\n" << Adense << std::endl;
  /*
   *	test of initA
   */
  MatrixXd Atest(6, 3), Atest2(6, 3);
  // clang-format off
  Atest << -1., 1., 0., 
           -1., 0., 1., 
           -1., 0., 0., 
            0., -1., 1.,
            0., -1., 0.,
            0., 0., -1.;
  Atest2 << -1., 1., 0., 
           -1., 0., 1., 
            0., -1., 1.,
           -1., 0., 0., 
            0., -1., 0.,
            0., 0., -1.;
  // clang-format on
  
  if (Adense.size() == Atest.size()) {
    // Test two possible row orderings, more options exist.
    if (Adense == Atest || Adense == Atest2)
      std::cout << "Test for initA passed!\n\n";
    else
      // NOTE: Other ordering of rows is also possible. 
      std::cout << "Test for initA failed: wrong output.\n"
                << "Note: only two possible row orderings tested, more options exist (see main.cpp).\n\n";
  } else
    std::cout << "Test for initA failed: wrong size.\n\n";

  /*
   *	run solveExtendedNormalEquations
   */
  MatrixXd D(4, 4);
  D << 0.0, -3.0, -4.0, -2.0, 3.0, 0.0, -1.0, 1.0, 4.0, 1.0, 0.0, 2.0, 2.0,
      -1.0, -2.0, 0.0;

  VectorXd x = solveExtendedNormalEquations(D);
  std::cout << "Solution extended normal equations:\n" << x << std::endl;
  /*
   *	test of solveExtendedNormalEquations
   */
  VectorXd xtest(3);
  xtest << 2, -1, -2;

  if (x.size() == xtest.size()) {
    if ((x - xtest).norm() < 1e-10)
      std::cout << "Test for solveExtendedNormalEquations passed!\n\n";
    else
      std::cout
          << "Test for solveExtendedNormalEquations failed: wrong output.\n\n";
  } else
    std::cout
        << "Test for solveExtendedNormalEquations failed: wrong size.\n\n";

  /*
   *	run solveNormalEquations
   */
  VectorXd y = solveNormalEquations(D);
  std::cout << "Solution normal equations:\n" << y << std::endl;
  /*
   *	test of solveNormalEquations
   */
  if (y.size() == xtest.size()) {
    if ((y - xtest).norm() < 1e-10)
      std::cout << "Test for solveNormalEquations passed!\n\n";
    else
      std::cout << "Test for solveNormalEquations failed: wrong output.\n\n";
  } else
    std::cout << "Test for solveNormalEquations failed: wrong size.\n\n";
  /*
   */ 
}
