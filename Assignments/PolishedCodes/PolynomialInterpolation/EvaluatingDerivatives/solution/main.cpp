#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#include "eval_deriv.hpp"

int main() {
  std::cout << "\nEnter \"0\" to test all functions." << std::endl
            << "Enter \"1\" to only test evaldp()." << std::endl
            << "Enter \"2\" to only test evaldp_naive()." << std::endl
            << "Enter \"3\" to only test polyTestTime()." << std::endl
            << "Enter \"4\" to only test dipoleval()." << std::endl
            << "Enter \"5\" to only test dipoleval_alt()." << std::endl
            << "Enter \"6\" to only test testDipolEval()." << std::endl
            << "Enter \"7\" to only test plotPolyInterpolant()." << std::endl;

  int ans = 0;
  std::cin >> ans;

  // Test evaldp()
  if (ans == 1 || ans == 0) {
    double x = -2.2;
    Eigen::VectorXd c(5);
    c << -1, -4, 2, 9, 4;
    // std::vector<double> c{-1,-4,2,9,4};
    std::pair<double, double> p = evaldp(c, x);
    std::cout << "Testing evaldp():" << std::endl;
    std::cout << "(p(x), p'(x)) = ( " << p.first << ", " << p.second << " )"
              << std::endl
              << std::endl;
  }

  // Test evaldp_naive()
  if (ans == 2 || ans == 0) {
    double x = -2.2;
    Eigen::VectorXd c(5);
    c << -1, -4, 2, 9, 4;

    std::pair<double, double> p_naive = evaldp_naive(c, x);

    std::cout << "Testing evaldp_naive():" << std::endl;
    std::cout << "(p(x), p'(x)) = ( " << p_naive.first << ", " << p_naive.second
              << " )" << std::endl
              << std::endl;
  }

  // Test polyTestTime()
  if (ans == 3 || ans == 0) {
    bool test1 = polyTestTime(15);
    std::cout << "polyTestTime() test passed: " << test1 << std::endl
              << std::endl;
  }

  // Test dipoleval()
  if (ans == 4 || ans == 0) {
    int n = 4;
    std::srand(199);
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(n, -1, 1);
    Eigen::VectorXd y = Eigen::VectorXd::Random(n);
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(5, -2, 2);

    std::cout << "Testing dipoleval():\n";
    std::cout << "p'(x) = " << dipoleval(t, y, x).transpose() << std::endl
              << std::endl;
  }

  // Test dipoleval_alt()
  if (ans == 5 || ans == 0) {
    int n = 4;
    std::srand(199);
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(n, -1, 1);
    Eigen::VectorXd y = Eigen::VectorXd::Random(n);
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(5, -2, 2);

    std::cout << "Testing dipoleval_alt():" << std::endl;
    std::cout << "p'(x) = " << dipoleval_alt(t, y, x).transpose() << std::endl
              << std::endl;
  }

  // Test testDipolEval()
  if (ans == 6 || ans == 0) {
    bool test2 = testDipolEval();
    std::cout << "testDipolEval() test passed: " << test2 << std::endl
              << std::endl;
  }

  // Test plotPolyInterpolant()
  if (ans == 7 || ans == 0) {
    plotPolyInterpolant("cos_plot");
    std::cout << "Running plotPolyInterpolant() ..." << std::endl
              << "Generated Plot." << std::endl;
  }

  return 0;
}
