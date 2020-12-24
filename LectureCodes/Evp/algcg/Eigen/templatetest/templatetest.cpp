
#include <Eigen/Dense> 
#include <iostream>

template <typename Function, typename T>  // Q: the T has to be the argument of Function??

void dummy(Function fu, T vec, double dum){

    unsigned int n = vec.size(); // specify in function that T has to be vector?
    for( unsigned int i = 0; i < n; ++i ){
        std::cout << vec[i] << std::endl ;}

}; 
