//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
////
#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "chebpolmult.hpp"

using namespace std;

// Compute the best approximation of the function $f$
// with Chebyshev polynomials.
// $\alpha$ is the output vector of coefficients.
template <typename Function>
void bestpolchebnodes(const Function &f, Eigen::VectorXd &alpha) {
    int n=alpha.size()-1;
    Eigen::VectorXd fn(n+1);
    for (int k=0; k<n+1; k++) {
        double temp=cos(M_PI*(2*k+1)/2/(n+1));
        fn(k)=f(temp);
    }
    
    vector<double> V;
    Eigen::MatrixXd scal(n+1,n+1);
    for (int j=0; j<n+1; j++) {
        V=chebpolmult(n,cos(M_PI*(2*j+1)/2/(n+1)));
        for (int k=0; k<n+1; k++) scal(j,k)=V[k];
    }
    
    for (int k=0; k<n+1; k++) {
        alpha(k)=0;
        for (int j=0; j<n+1; j++) {
            alpha(k)+=2*fn(j)*scal(j,k)/(n+1);
        }
    }
    alpha(0)=alpha(0)/2;
}

int main(){
	int n;
	
// Check the orthogonality of Chebyshev polynomials
    n=10;
    vector<double> V;
    Eigen::MatrixXd scal(n+1,n+1);
    for (int j=0; j<n+1; j++) {
        V=chebpolmult(n,cos(M_PI*(2*j+1)/2/(n+1)));
        for (int k=0; k<n+1; k++) scal(j,k)=V[k];
    }
    cout<<"Scalar products: "<<endl;
    for (int k=0; k<n+1; k++)
        for (int l=k+1; l<n+1; l++)
            cout<<scal.col(k).dot(scal.col(l))<<endl;

// Test the implementation
    auto f = [] (double & x) {return 1/(pow(5*x,2)+1);};
    n=20;
    Eigen::VectorXd alpha(n+1);
    bestpolchebnodes(f, alpha);
    
    // Compute the error
    Eigen::VectorXd X = Eigen::VectorXd::LinSpaced(1e6,-1,1);
    auto qn = [&alpha,&n] (double & x) {
        double temp=0;
        vector<double> V=chebpolmult(n,x);
        for (int k=0; k<n+1; k++) temp+=alpha(k)*V[k];
        return temp;
    };
    double err_max=abs(f(X(0))-qn(X(0)));
    for (int i=1; i<1e6; i++) err_max=std::max(err_max,abs(f(X(i))-qn(X(i))));
    cout<<"Error: "<< err_max <<endl;
}
