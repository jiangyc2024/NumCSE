#include<Eigen/Dense>
#include "laserquad.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp; 

int main() {
    auto f = [] (double x) { return x*x*x*std::log(std::abs(x)+1); };
    QuadRule Q;
    Q.nodes.resize(5);
    Q.nodes << -1., -std::sqrt(3./7.), 0, std::sqrt(3./7.), 1.;
    Q.weights.resize(5);
    Q.weights << 0.1, 49./90., 32./45., 49./90., 0.1;
    std::cout << "Test of evalquad() for log(|x|+1)*x^3: "
              << evalquad(-2,3,f,Q) << "\n\n";
    
    
    auto g = [] (double x, double y) { return x*(y-1.)*std::sin(6.*x*y); };
    std::cout << "Test of gaussquadtriangle() for x*(y-1)*sin(6xy): "
              << gaussquadtriangle(g,5) << "\n\n";
    
    convtest2DQuad();
    return 0;
}
