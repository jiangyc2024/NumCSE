#include "SDIRK.hpp"


int main() {
  double rate = cvgSDIRK();
  std::cout << std::endl << "The rate is " << rate << std::endl;
}
