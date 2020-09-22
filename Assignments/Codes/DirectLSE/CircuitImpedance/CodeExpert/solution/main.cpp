#include <iostream>

#include "circuitimpedance.hpp"

int main() {
	NodalPotentials(1, 1);
	ImpedanceMap test(1, 1);
	std::cout << test(1) << " " << test(2) << " " << test(2.5) << std::endl;
}
