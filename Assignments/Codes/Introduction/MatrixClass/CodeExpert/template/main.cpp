
#include "MatrixClass.hpp"

int main(){
    
    std::cout<< "\nTest of smallTriangular():\n";
    std::cout<< smallTriangular(1,2,3) << std::endl;
    std::cout<< "\nTest of constantTriangular():\n";
    std::cout<< constantTriangular(3,20) << std::endl;
    std::cout<< "\nTest of arithmetics():\n";
    std::srand(5); // So that random behaviour is predictable.
    std::cout<< arithmetics(2,5) << std::endl;
    std::cout<< "\nTest of casting():\n";
    std::cout<< casting() << std::endl;    

    return 0;
    }
