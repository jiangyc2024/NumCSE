
#include "tridiagqr.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

using namespace Eigen;

int main() {
    int n = 6;
    VectorXd d(n), l(n-1), u(n-1);
    l = VectorXd::LinSpaced(n-1,1,n-1);
    d = VectorXd::LinSpaced(n,n,2*n-1);
    u = VectorXd::LinSpaced(n-1,2*n,3*n-2);
    MatrixXd Adense(n,n);
    Adense.setZero();
    Adense.diagonal() = d;
    Adense.diagonal(-1) = l;
    Adense.diagonal(1) = u;
    std::cout << std::setprecision(3);
    std::cout << "A=\n" << Adense << std::endl;
    
    TriDiagonalMatrix A = TriDiagonalMatrix(d,l,u);
    TriDiagonalQR Aqr = TriDiagonalQR(A);
    
    MatrixXd Q,R;
    Aqr.getDense(Q,R);
    std::cout << "\nQ=\n" << Q << std::endl;
    std::cout << "R=\n" << R << std::endl;
    
    std::cout << "Q orthonormal, |Q^TQ-Id|=\n" << (Q.transpose()*Q-MatrixXd::Identity(n,n)).norm()<< std::endl;
    std::cout << "|QR-A|=\n" << (Q*R - Adense).norm() << std::endl;
    std::cout << "|Q^TA-R|=\n" << (Q.transpose()*Adense - R).norm() << std::endl;
    
    VectorXd x,y;
    x = VectorXd::Random(n);
    VectorXd QTxdense = Q.transpose()*x;
    y = Aqr.applyQT(x);    
    std::cout << "|dense - tridiag|=\n" << (QTxdense-y).norm() << std::endl;
    
    // Solve Ax=b
    VectorXd b = VectorXd::Random(n);
    x = Aqr.solve(b);
    y = Adense.lu().solve(b);
    std::cout << "|dense - tridiag|:\n" << (x-y).norm() << std::endl;
    
    
    if(false) {
        std::cout << "\n Givens params\n";
        Vector2d a(2);
        double rho, gamma, sigma;
        std::tuple<double,double,double> params;
        Matrix2d G, Grho;
        
        std::cout << "\n Givens matrix and rotated vector\n";
        a << -1,1;
        params = compGivensRotation(a);
        rho = std::get<0>(params);
        gamma = std::get<1>(params);
        sigma = std::get<2>(params);
        std::cout << "rho = " << rho << std::endl;
        G << gamma, sigma, -sigma, gamma;
        Grho = Givens(rho);
        std::cout << G << " *\n" << a << " =\n" << G.transpose()*a << std::endl;
        std::cout << "|G-Grho|: " << (G-Grho).norm() << std::endl;
        
        
        std::cout << "\n Givens matrix and rotated vector\n";
        a << 0,2;
        params = compGivensRotation(a);
        rho = std::get<0>(params);
        gamma = std::get<1>(params);
        sigma = std::get<2>(params);
        std::cout << "rho = " << rho << std::endl;
        G << gamma, sigma, -sigma, gamma;
        Grho = Givens(rho);
        std::cout << G << " *\n" << a << " =\n" << G.transpose()*a << std::endl;
        std::cout << "|G-Grho|: " << (G-Grho).norm() << std::endl;
        
        std::cout << "\n Givens matrix and rotated vector\n";
        a << 2,0;
        params = compGivensRotation(a);
        rho = std::get<0>(params);
        gamma = std::get<1>(params);
        sigma = std::get<2>(params);
        std::cout << "rho = " << rho << std::endl;
        G << gamma, sigma, -sigma, gamma;
        Grho = Givens(rho);
        std::cout << G << " *\n" << a << " =\n" << G.transpose()*a << std::endl;
        std::cout << "|G-Grho|: " << (G-Grho).norm() << std::endl;
        
        std::cout << "\n Givens matrix and rotated vector\n";
        a << 4,3;
        params = compGivensRotation(a);
        rho = std::get<0>(params);
        gamma = std::get<1>(params);
        sigma = std::get<2>(params);
        std::cout << "rho = " << rho << std::endl;
        G << gamma, sigma, -sigma, gamma;
        Grho = Givens(rho);
        std::cout << G << " *\n" << a << " =\n" << G.transpose()*a << std::endl;
        std::cout << "|G-Grho|: " << (G-Grho).norm() << std::endl;
    }
    return 0;
}
