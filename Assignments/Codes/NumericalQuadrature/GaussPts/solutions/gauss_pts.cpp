#include <vector>
#include <iostream>
#include <Eigen/Dense>

//! \brief Golub-Welsh implementation 5.3.35
//! \param[in] n number of Gauss nodes
//! \param[out] w weights for interval [-1,1]
//! \param[out] xi ordered nodes for interval [-1,1]
void gaussrule(int n, Eigen::VectorXd & w, Eigen::VectorXd & xi) {
    w.resize(n);
    xi.resize(n);
    if( n == 0 ) {
        xi(0) = 0;
        w(0) = 2;
    } else {
        Eigen::VectorXd b(n-1);
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n,n);
        
        for(int i = 1; i < n; ++i) {
            double d = (i) / sqrt(4. * i * i - 1.);
            J(i,i-1) = d;
            J(i-1,i) = d;
        }
        
        Eigen::EigenSolver<Eigen::MatrixXd> eig(J);
        
        xi = eig.eigenvalues().real();
        w = 2 * eig.eigenEigen::VectorXds().real()
                                  .topRows<1>()
                                  .cwiseProduct( 
                                      eig.eigenEigen::VectorXds()
                                         .real()
                                         .topRows<1>()
                                  );
    }
    
    std::Eigen::VectorXd<std::pair<double,double>> P;
    P.reserve(n);
    for(unsigned int i = 0; i < n; ++i) {
        P.push_back(std::pair<double,double>(xi(i),w(i)));
    }
    std::sort(P.begin(), P.end());
    for(unsigned int i = 0; i < n; ++i) {
        xi(i) = std::get<0>(P[i]);
        w(i) = std::get<1>(P[i]);
    }
}

//! \brief Compute the function g in the Gauss nodes
//! \param[in] f object with an evaluation operator (e.g. a lambda function) representing the function f
//! \param[in] n number of nodes
//! \param[out] Eigen::VectorXd containing the function g calculated in the Gauss nodes.
/* SAM_LISTING_BEGIN_1 */
template<typename Function>
Eigen::VectorXd comp_g_gausspts(Function f, unsigned int n) {

    Eigen::VectorXd g(n);
    Eigen::VectorXd w(n), xi(n);
    gaussrule(n, w, xi);  // Compute Gauss nodes and weights relative to [-1,1]

    Eigen::VectorXd p(n), q(n), expxi(n), fxi(n);
    w = w/2.;  // Rescale the weights to [0,1]
    xi = (xi + Eigen::MatrixXd::Ones(n,1))/2.;  // Rescale the nodes to [0,1]
    
    for (int l=0; l<n; l++) {
        fxi(l) = f(xi(l));
        expxi(l) = exp(xi(l));
    }
    
    // Initialize the Eigen::VectorXds p and q
    p(0) = 0;
    q(n-1) = w(n-1) * expxi(n-1) * fxi(n-1);
    
    // Fill in the Eigen::VectorXds p and q (O(n) complexity)
    for (int l=0; l<n-1; l++) {
        p(l+1) = p(l) + w(l)/expxi(l)*fxi(l);
        q(n-2-l) = q(n-1-l) + w(n-2-l)*expxi(n-2-l)*fxi(n-2-l);
    }
    
    // Finally, constructing the output (O(n) complexity)
    g = (p.array() * expxi.array() + q.array() / expxi.array()).matrix();
    return g;
}
/* SAM_LISTING_END_1 */


int main() {
    /* SAM_LISTING_BEGIN_2 */
    int n = 21;

    auto f =  [] (double y) {return std::exp(-std::abs(.5-y)); };

    Eigen::VectorXd g = comp_g_gausspts(f,n);
    std::cout << "g(xi^" << n << "_" << (n+1)/2 << ") = " << g((n+1)/2-1) << std::endl;
    /* SAM_LISTING_END_2 */
}
