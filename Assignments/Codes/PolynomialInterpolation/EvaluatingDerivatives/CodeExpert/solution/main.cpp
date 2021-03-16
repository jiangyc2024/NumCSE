
#include <iomanip>
#include <iostream>
#include <vector>

#include "eval_deriv.hpp"

int main() {
  std::cout << "\nEnter \"0\" to test all functions.\n"
            << "Enter \"1\" to only test evaldp().\n"
            << "Enter \"2\" to only test evaldp_naive().\n"
            << "Enter \"3\" to only test polyTestTime().\n"
            << "Enter \"4\" to only test dipoleval().\n"
            << "Enter \"5\" to only test dipoleval_alt().\n"
            << "Enter \"6\" to only test testDipolEval().\n"
            << "Enter \"7\" to only test plotPolyInterpolant().\n\n";
  
  int ans=0;
  std::cin >> ans;
  
  if (ans==1 || ans==0) {
    double x = -2.2;
    std::vector<double> c{-1,-4,2,9,4};
    std::pair<double,double> p = evaldp(c,x);
    std::cout << "Testing evaldp():\n";
    std::cout << "(p(x), p'(x)) = ( " << p.first << ", " << p.second << " )\n\n";
  }
  if (ans==2 || ans==0) {
    double x = -2.2;
    std::vector<double> c{-1,-4,2,9,4};
    std::pair<double,double> p_naive = evaldp_naive(c,x);
    std::cout << "Testing evaldp_naive():\n";
    std::cout << "(p(x), p'(x)) = ( " << p_naive.first << ", " << p_naive.second << " )\n\n";
  }
  if (ans==3 || ans==0) {
    bool test1 = polyTestTime(15);
    std::cout << "polyTestTime() test passed: " << test1 << "\n\n";
  }
  if (ans==4 || ans==0) {
    int n = 4;
    std::srand(199);
    VectorXd t = VectorXd::LinSpaced(n,-1,1);
    VectorXd y = VectorXd::Random(n);
    VectorXd x = VectorXd::LinSpaced(5,-2,2);
    std::cout << "Testing dipoleval():\n";
    std::cout << "p'(x) = " << dipoleval(t, y, x).transpose() << "\n\n";
  }
  if (ans==5 || ans==0) {
    int n = 4;
    std::srand(199);
    VectorXd t = VectorXd::LinSpaced(n,-1,1);
    VectorXd y = VectorXd::Random(n);
    VectorXd x = VectorXd::LinSpaced(5,-2,2);
    std::cout << "Testing dipoleval_alt():\n";
    std::cout << "p'(x) = " << dipoleval_alt(t, y, x).transpose() << "\n\n";
  }
  if (ans==6 || ans==0) {
    bool test2 = testDipolEval();
    std::cout << "testDipolEval() test passed: " << test2 << "\n\n";
  }
  if (ans==7 || ans==0) {
    plotPolyInterpolant("cos_plot");
  }
  
  return 0;
}
