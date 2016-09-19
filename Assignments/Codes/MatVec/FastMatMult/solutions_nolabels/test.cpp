//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>

#include "strassen.cpp"

int main()
{
    // seed random number generator
    srand((unsigned int) time(0));
    
    // check algorithm for correctness
    int k=2;
    int n=pow(2,k);
    MatrixXd A=MatrixXd::Random(n,n);
    MatrixXd B=MatrixXd::Random(n,n);
    MatrixXd AB(n,n), AxB(n,n);
    AB=strassenMatMult(A,B);
    AxB=A*B;
    std::cout << "Using Strassen's method, A*B=" << AB << std::endl;
    std::cout << "Using standard method, A*B=" << AxB << std::endl;
    std::cout << "The norm of the error is " << (AB-AxB).norm() << std::endl;
}
