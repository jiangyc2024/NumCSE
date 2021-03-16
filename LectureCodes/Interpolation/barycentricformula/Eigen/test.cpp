///////////////////////////////////////////////////////////////////////////
/// Testing code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Ralf Hiptmair
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#include "intpolyval.hpp"
#include "ipvclass.hpp"
#include <iostream>
#include <list>

using namespace std;
using namespace Eigen;

int main() {
  int n = 5;
  int N = 7;
  VectorXd t = VectorXd::LinSpaced(n, 1, n);
  VectorXd x = VectorXd::LinSpaced(N, 0.0, 1.0);
  VectorXd y = t.cwiseSqrt();
  VectorXd p(N);


  std::cout << "Polynomial interpolation" << std::endl;
  std::cout << "Nodes = " << t.transpose() << std::endl;
  std::cout << "Values = " << y.transpose() << std::endl;
  std::cout << "Evaluation points = " << x.transpose() << std::endl;
  
  // Evaluation by means of function intpolyval()
  intpolyval(t, y, x, p);

  // Evaluation by means of interpolation class
  BarycPolyInterp<> ipc(t);
  cout << "p = " << p.transpose() << endl
       << "ipc.p = " << ipc.eval<VectorXd>(y, x).transpose() << endl;

  // Construction from STL
  std::list<std::complex<double>> lst;
  for (int i = 0; i < n; i++)
    lst.push_back(t(i));
  BarycPolyInterp<std::complex<double>> ipc2(lst);
  VectorXcd xc = x.cast<complex<double>>();
  cout << "ipc.p = " << ipc2.eval<VectorXcd>(y, xc).transpose() << endl;

  return 0;
}
