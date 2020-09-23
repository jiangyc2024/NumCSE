#include <iomanip>
#include <iostream>

#include "circuitimpedance.hpp"

int main() {
  // Create a NodalPotentials object with R = 1, Rx = 1
  NodalPotentials NP(1, 1);
  // Print NodalPotentials for some source voltages
  std::cout << "--> Some NodalPotentials" << std::endl;
  std::cout << "Nodal potentials for V = 1:" << std::endl;
  std::cout << NP(1) << std::endl << std::endl;
  std::cout << "Nodal potentials for V = 2:" << std::endl;
  std::cout << NP(2) << std::endl << std::endl;
  std::cout << "Nodal potentials for V = 2.5:" << std::endl;
  std::cout << NP(2.5) << std::endl << std::endl;

  // Create an ImpedanceMap with resistance 1 and voltage 1
  ImpedanceMap IM(1, 1);

  // Print a table with various impedances
  std::cout << "--> Table of impedances" << std::endl;
  // Table header
  std::cout << std::setw(30) << "Impedance [Ohm]" << std::setw(30)
            << "R_x [Ohm]" << std::endl;
  // Table content: print impedance for various resistance values
  for (auto Rx = 1; Rx <= 1024; Rx *= 2) {
    std::cout << std::setw(30) << IM(Rx) << std::setw(30) << Rx << std::endl;
  }
}
