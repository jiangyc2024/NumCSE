/***********************************************************************
 *                                                                     *
 * Demo code                                                           *    
 * (Prof. Dr. R. Hiptmair)                                             * 
 * Author: R.H.                                                        *
 * Date: May 2018                                                      * 
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

// Header for basic IO
#include <iostream>
// Provides random acccess container class
#include <vector>
// Provides function types 
#include <functional>

// Definition of functions
double binop(double arg1,double arg2) { return (arg1/arg2); }

void stdfunctiontest() {
  // Vector of objects of a particular signature
  std::vector<std::function<double(double,double)>> fnvec;
  // Store reference to a regular function
  fnvec.emplace_back(binop);
  // Store are lambda function
  fnvec.emplace_back([] (double x,double y)  -> double { return y/x; });
  for (auto const & fn : fnvec) { std::cout << fn(3,2) << std::endl; }
}

int main() {
  // Call testing function;
  stdfunctiontest();
  return (0);
}

