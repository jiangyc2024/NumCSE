#include <iostream>
#include <iomanip>

#include "circuitimpedance.hpp"

int main() {
	NodalPotentials(1, 1);
	
	// Create a first ImpedanceMap with resistance 1 and voltage 1
	ImpedanceMap IM = ImpedanceMap(1, 1);
	// Test Impedance at value one
	std::cout << "--> Impedance with R = 1: " << IM(1) << std::endl;
	
	// Print a table with various impedances
	std::cout << "--> Table of impedances" << std::endl;
	// Table header
	std::cout << std::setw(30) << "Impedance [Ohm]"
	<< std::setw(30) << "R_x [Ohm]"
	<< std::endl;
	// Table content: print impedance for various resistance values
	std::cout << std::setw(30) << IM(0.1)
	<< std::setw(30) << 0.1
	<< std::endl;
	for(auto Rx = 1; Rx <= 1024; Rx *= 2) {
		std::cout << std::setw(30) << IM(Rx)
		<< std::setw(30) << Rx
		<< std::endl;
	}
}
