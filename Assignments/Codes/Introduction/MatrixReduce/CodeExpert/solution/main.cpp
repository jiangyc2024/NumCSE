

#include "MatrixReduce.hpp"
#include <iostream>


int main(){
    
    std::cout<< "\nTest of average():\n";
    Matrix3d A = Matrix3d::Identity();
    std::cout<< average( A ) << "\n\n";
    std::cout<< "Test of percent_zero():\n";
    std::cout<< percent_zero( A ) << "\n\n";
    
    std::cout<< "Test of has_zero_column():\n";
    MatrixXd B = MatrixXd::Identity(4,5);
    std::cout<< has_zero_column( B ) << " " << has_zero_column( B.transpose() ) << "\n\n";
    
    std::cout<< "Test of columns_sum_to_zero():\n";
    std::srand(5); // So that random behaviour is predictable.
    Matrix3d C = Matrix3d::Random() + Matrix3d::Constant(1);
    std::cout<< columns_sum_to_zero( C ) << "\n\n";
    
    return 0;
    }
