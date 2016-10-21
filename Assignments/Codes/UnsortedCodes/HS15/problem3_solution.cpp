#include <iostream>

#include <vector>

#include <Eigen/Dense>

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

//! \brief Function for the evaluation of minmod(x,y)
double minmod(double x, double y) {
    //// TODO: problem 3a, implement this function
    
    return x*y <= 0. ? 0. : (x >= 0. ? std::min(x,y) : std::max(x,y)); 
}

//! \brief Function for the computation of the values \f$ \Delta_j \f$
//! \param[in] t Vector of nodes \f$ t_j, j = 0,\dots,n \f$ 
//! \param[in] y Vector of values \f$ y_j, j = 0,\dots,n \f$ 
//! \return Vector \f$ (\Delta_1, \dots, \Delta_n) \f$ relative to (t, y) data
Vector delta(const Vector & t, const Vector & y) {
    assert(t.size() == y.size());
    assert(t.size() > 1);
    
    //// TODO: problem 3b, implement this function
    
    Vector dlt(t.size()-1);
    
    for(unsigned int j = 0; j < t.size()-1; ++j) {
        dlt(j) = (y(j+1) - y(j)) / (t(j+1) - t(j));
    }
     
    return dlt;
}

//! \brief Class representing a pecevise cubic Hermite interpolant \f$ s \f$ with minmod-reconstructed slopes
class mmPCHI {
public:
    //! \brief Construct the slopes \f$ c_j \f$ of the interpolant, from the data \f$ (t_j, y_j) \f$.
    //! \param[in] t Vector of nodes \f$ t_j, j = 0,\dots,n \f$ 
    //! \param[in] y Vector of values \f$ y_j, j = 0,\dots,n \f$ 
    mmPCHI(const Vector &t, const Vector &y);
    
    //! \brief Evaluate interpolant defined by (*this) at the points defined by x and store the result in v
    //! \param[in] x Vector of points \f$ x_i, i = 0,\dots,m \f$ 
    //! \return v Vector of values \f$ s(x_i), i = 0,\dots,m \f$ 
    Vector eval(const Vector &x) const;
private:
    const Vector t, y;//!< Vectors containing nodes and values
    Vector c; //!< Vector containing slopes for the interpolant
};

mmPCHI::mmPCHI(const Vector &t, const Vector &y)
: t(t), y(y), c(t.size()) {
    
    assert(t.size() == y.size());
    assert(t.size() > 1);
    
    //// TODO: problem 3c, implement this function
    
    unsigned int n = t.size()-1;
    
    Vector D = delta(t,y);
    
    // Reconstruct slopes using minmod
    c(0) = D(0);
    for(unsigned int i = 1; i < n; ++i) {
        
        c(i) = minmod(D(i), D(i-1));
    }
    c(n) = D(n-1);
}

Vector mmPCHI::eval( const Vector &x) const {
    Vector v(x.size());
    
    //// TODO: problem 3d, implement this function
    Vector xcopy = x;
    std::sort(xcopy.data(), xcopy.data()+xcopy.size());
    
    unsigned int i = 0;
    for(unsigned int j = 0; j < t.size()-1; ++j) {
        double t1 = t(j);
        double t2 = t(j+1);
        double y1 = y(j); 
        
        double h = t2 - t1;
        double a1 = y(j+1) - y1;
        double a2 = a1 - h*c(j);
        double a3 = h*c(j+1) - a1 - a2;
        while( i < xcopy.size() && xcopy(i) <= t2 ) {
            double tx = (xcopy(i) - t1) / h;
            v(i) = y1 + (a1+(a2+a3*tx)*(tx-1.))*tx;
            i++;
        }   
    }
    
    return v;
}

int main(int, char**) {
    //// PROBLEM 3 TEST
    
    auto f = [] (double t) -> double { return 1. / (1. + t*t); };
    //auto f = [] (double t) -> double { return 1. / (1. + std::abs(t)); };
    
    // Interval
    const double a = -5., b = 5;
    // Number of sampling nodes
    const unsigned int M = 10000;
    
    // Sampling nodes and values
    const Vector x = Vector::LinSpaced(M,a,b);
    Vector v_ex(M);
    // Fill in exact values of f
    for(unsigned int i = 0; i<M; ++i) {
        v_ex(i) = f(x(i));
    }
    // Output table with error w.r.t number of nodes
    std::cout << "n" << "\t" << "L^inf err." << std::endl;
    for(unsigned int n = 4; n <= 2048; n=n<<1) {
        Vector t = Vector::LinSpaced(n,a,b);
        Vector y(n);
        for(unsigned int i = 0; i<n; ++i) {
            y(i) = f(t(i));
        }
        
        mmPCHI s(t,y);
        const Vector v = s.eval(x);
        std::cout << n << "\t\t" << (v - v_ex).lpNorm<Eigen::Infinity>() << std::endl;
    }
    
}
 
