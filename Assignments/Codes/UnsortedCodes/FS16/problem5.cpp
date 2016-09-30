#include <iostream>

#include <vector>

#include <Eigen/Dense>

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

//! \brief Evolves the vector y0 using the evolution operator \tilde{\Psi} with step-size y0
//! \tparam DiscEvlOp type for evolution operator (e.g. lambda function type)
//! \param[in] Psi original evolution operator, must have operator(double, const Vector&)
//! \param[in] p parameter p for construction of Psi tilde
//! \param[in] h step-size
//! \param[in] y0 previous step
//! \return Evolved step \f$ \tilde{\Psi}^h y0 \f$
template <class DiscEvlOp>
Vector psitilde(const DiscEvlOp& Psi, unsigned int p, double h, const Vector & y0) {
    // TODO: problem 5b: implement psi tilde operator
}

//! \brief Evolves the vector y0 using the evolution operator from time 0 to T using equidistant steps
//! \tparam DiscEvlOp type for evolution operator (e.g. lambda function type)
//! \param[in] Psi original evolution operator, must have operator(double, const Vector&)
//! \param[in] T final time
//! \param[in] y0 initial data
//! \param[in] N number of steps
//! \return Vector of all steps y_0, y_1, ...
template <class DiscEvlOp>
std::vector<Vector> odeintequi(const DiscEvlOp& Psi, double T, const Vector &y0, int N) {
    
    std::vector<Vector> Y;
    
    // TODO: problem 5c: implement equidistant solver for evolution operaotr Psi
    
    return Y;
}

//! \brief Evolves the vector y0 using the evolution operator from time 0 to T using adaptive error control
//! \tparam DiscEvlOp type for evolution operator (e.g. lambda function type)
//! \param[in] Psi original evolution operator, must have operator(double, const Vector&)
//! \param[in] T final time
//! \param[in] y0 initial data
//! \param[in] h0 initial step size
//! \param[in] p parameter p for construction of Psi tilde
//! \param[in] reltol relative tolerance for error control
//! \param[in] abstol absolute tolerance for error control
//! \param[in] p parameter p for construction of Psi tilde
//! \param[in] hmin minimal step size
//! \return Vector of all steps y_0, y_1, ...
template <class DiscEvlOp>
std::vector<Vector> odeintssctrl(const DiscEvlOp& Psi, double T, const Vector &y0, double h0,
                                 unsigned int p, double reltol, double abstol, double hmin) {
    double t = 0.;
    double h = h0;
    std::vector<Vector> Y;
    Y.push_back(y0);
    Vector y = y0;
    
    while( t < T && h > hmin ) {
        Vector y_high = Vector::Zero(y0.size()); // TODO: problem 5e: fix this line
        Vector y_low = Vector::Zero(y0.size()); // TODO: problem 5e: fix this line
        double err_est = (y_high - y_low).norm();
        
        if( true /* TODO: problem 5e: fix this line */ ) {
            y = y_high;
            Y.push_back(y_high);
            
            t += std::min(T-t,h);
            h *= 1.1;
        } else {
            h /= 2.;
        }
    }

    return Y;
}

int main(int, char**) {
    
    auto f = [] (const Vector &y) -> Vector { return Vector::Ones(1) + y*y; };
    Vector y0 = Vector::Zero(1);
    double T = 1.;
    auto y_ex = [] (double t) -> Vector { Vector y(1); y << tan(t); return y; };
    
    unsigned p = 1;
    
    //// PROBLEM 5d
    std::cout << "*** PROBLEM 5d:" << std::endl;
    
    // TODO: problem 5d: determine order of convergence of PsiTilde using odeintequi
    
    //// PROBLEM 5e
    std::cout << "*** PROBLEM 5e:" << std::endl;
    
    double h0 = 1./100.;
    std::vector<Vector> Y;
    double err = 0;
    
    // TODO: problem 5e: uncomment following lines to test odeintssctrl
//     Y = odeintssctrl(Psi, T, y0, h0, p, 10e-4, 10e-4, 10e-5);
//     err = (Y.back() - y_ex(1)).norm();
    std::cout << "Adaptive error control results:" << std::endl;
    std::cout << "Error" << "\t\t" << "No. Steps" << "\t" << "y(1)" << "\t" << "y(ex)" << std::endl;
//     std::cout << err << "\t" << Y.size() << "\t" << Y.back() << "\t" << y_ex(1) << std::endl;
}
                                
