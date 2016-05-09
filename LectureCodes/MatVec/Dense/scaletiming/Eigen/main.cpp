#define NDEBUG true
#include <Eigen/Dense>
#include <figure/figure.hpp>
#include <iostream>
#include <limits>
#include <cmath>
#include "timer.hpp"
using namespace std::chrono;
using namespace Eigen;
#include "scaletiming.hpp"
int main(){
  scaletiming();
  return 0;
}
