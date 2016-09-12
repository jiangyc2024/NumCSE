#include <iostream>
#include <Eigen/Dense>

int main () {
	std::cout << gausstiming() << std::endl;
	return 0;
}

/*Compiling options:
 * g++ -std=c++11  -m64 -I${MKLROOT}/include main.cpp  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl

 * variable MKLROOT must be set!

 */
