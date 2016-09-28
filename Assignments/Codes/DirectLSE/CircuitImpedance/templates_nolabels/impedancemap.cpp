//// 
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch> 
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.
//// Report issues to: https://gitlab.math.ethz.ch/NumCSE/NumCSE/issues
////
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/LU>


using namespace Eigen;

// Models a resistance between node i,j with i < j
// (in principle, it can handle different resistances)
typedef std::pair< int, int>            resistor;
// Vector containing all resistances in
// the circuit (excluded the variable one, which is modelled afterwards)
typedef std::vector< resistor >         resistor_topology;
// Models a ground or source of voltage, indexes store the node
// connected to the ground/source through a resistance
typedef std::tuple<int, int, double>    voltage;
// Vector containing a voltage object for each node connected to sink or source.
typedef std::vector< voltage >          voltage_topology;

//! \brief Class implementing the topology of the circuit (cf. Figure)
//! Computes impedance of the entire circuit (between node 16 and 17) exploiting
//! the SMW formula for the inversion of low rank perturbations of a matrix A0,
//! whose factorization in known in advance
class ImpedanceMap {
    std::size_t nnodes; //< Number of nodes in the circuit
public:
    //! \brief Constructor: Build system matrix and r.h.s. and perform a LU decomposition
    //! The LU decomposition is stored in 'lu' and can be reused in the SMW formula
    //! to avoid expensive matrix solves for repeated usages of the operator()
    //! \param R Resistance (in Ohm) value of $R$
    //! \param W Source voltage $W$ at node 16 (in Volt), ground is set to 0V at node 17
    ImpedanceMap(double R, double W) : R(R), W(W) {
        // TODO: build the matrix A_0 and compute its lu factorization. Compute the
        // right hand side and store it in rhs.
    };

    //! \brief Compute the impedance given the resistance $R_x$
    //! Use SMW formula for low rank perturbations to reuse LU
    //! factorization.
    //! \param Rx Resistence $R_x > 0$ of the varistor between node 14 and 15
    //! \return impedance $W / I$ of the system $A_{R_x}$
    double operator()(double Rx) {
    // TODO: use SMW formula to compute the solution of A_{R_x} x = rhs
    // and compute impedance of the system
};

int main(void) {
    // Create a first ImpedanceMap with resistance 1 and voltage 1
    ImpedanceMap IM = ImpedanceMap(1, 1);
    // Test Impedance at value one
    std::cout << "--> Impedance with R = 1: " << IM(1) << std::endl;

    // Print a table with various impedances
    std::cout << "--> Table of impedances" << std::endl;
    // Table header
    std::cout << std::setw(30) << "Impedance [Ohm]"
              << std::setw(30) << "R_x [Ohm]"
              << std::endl;
    // Table content: print impedance for various resistance values
    std::cout << std::setw(30) << IM(0)
              << std::setw(30) << 0
              << std::endl;
    std::cout << std::setw(30) << IM(0.1)
              << std::setw(30) << 0.1
              << std::endl;
    for(auto Rx = 1; Rx <= 1024; Rx *= 2) {
        std::cout << std::setw(30) << IM(Rx)
                  << std::setw(30) << Rx
                  << std::endl;
    }


}
