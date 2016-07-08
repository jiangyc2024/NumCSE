#include <Eigen/Dense>
#include <Eigen/LU>

#include <iostream>
#include <iomanip>

#include "timer.h"

// Models a resistence between node i,j with i < j (could, in principle, handle different resistences)
typedef std::pair< int, int>            resistor;
// Vector containing all resistences in the circuit (excluded the variable one, which is modelled afterwards)
typedef std::vector< resistor >         resistor_topology;
// Models a ground or source of voltage, indexes store the node connected to the ground/source through a resistence
typedef std::tuple<int, int, double>    voltage;
// Vector containing a voltage object for each node connected to sink or source.
typedef std::vector< voltage >          voltage_topology;

//! \brief Class implementing the topology of the circuit (cf. Figure)
//! Computes impedance of the entire circuit (between node 16 and 17) exploiting
//! the SMW formula for the inversion of low rank perturbations of a matrix A0,
//! whose factorization in known in advance
class ImpedanceMap {
    static const std::size_t nnodes = 15; //< Compile time hard-coded size of circuit
public:
    //! \brief Constructor: Builds system matrix and rhs and performs lu decomposition
    //! The lu decomposition is stored in lu and can be reused in the SMW formula
    //! to avoid expensive matrix solves for repeated usages of the operator()
    //! \param R_ resistence R = R1 = R2 = R3 = R4 stored in R
    //! \param U_ source voltage U at node 16, gound is 0 at node 17
    ImpedanceMap(double R_, double U_) : R(R_), U(U_) {
        // In the following, instead of specifying directly the entries of A\_0 and of rhs
        // we define some auxiliary structure and automatically generate the right entries for
        // the structures we specified.
        // This way: we avoid silly mistakes, avoid a very long list of matrix entreis specifications,
        // have a flexible way (we may easily change the topology) and automatically take care of shifting the indices
        // Of course, if you want, you can manually put each entry in the matrix and rhs. The result (should) be the same)
        
        // We implement the topology of the resistances by pushing each
        // Resistance between node i < j in a std::vector of pair (i,j) \in \mathbb{N}^2
        // This way we can automatically build a symmetric matrix, take care of indexing 
        // and avoid dumb mistakes during the matrix filling
        resistor_topology T;
        T.reserve(24);
        T.push_back(resistor(1,2));
        T.push_back(resistor(1,5));
        T.push_back(resistor(2,5));
        T.push_back(resistor(2,3));
        T.push_back(resistor(2,14));
        T.push_back(resistor(3,4));
        T.push_back(resistor(3,15));
        T.push_back(resistor(4,6));
        T.push_back(resistor(4,15));
        T.push_back(resistor(5,7));
        T.push_back(resistor(5,14));
        T.push_back(resistor(6,9));
        T.push_back(resistor(6,15));
        T.push_back(resistor(7,10));
        T.push_back(resistor(7,11));
        T.push_back(resistor(8,9));
        T.push_back(resistor(8,12));
        T.push_back(resistor(8,13));
        T.push_back(resistor(8,15));
        T.push_back(resistor(9,13));
        T.push_back(resistor(10,11));
        T.push_back(resistor(11,12));
        T.push_back(resistor(12,13));
        // Default matrix entries sert to one (varistor)
        T.push_back(resistor(14,15));
        
        // We implement voltage gound and source topology by specifiying which node is connected to ground/source
        // trough a resistance. This will be also part of rhs
        voltage_topology S;
        S.reserve(4);
        S.push_back(voltage(6, 16, U));
        S.push_back(voltage(7, 17, 0));
        S.push_back(voltage(11,17, 0));
        S.push_back(voltage(14,17, 0));
        
        // Automatically build A0 filling from topology
        Eigen::Matrix<double, nnodes,nnodes> A0 = Eigen::Matrix<double,nnodes,nnodes>::Zero();
        for(auto it: T) {
            // Shift indeces down by 1
            auto i = it.first - 1;
            auto j = it.second - 1;
            // Fill diagonal: each resistance contributes += R into the diagonal of both nodes
            A0(i,i) += 1;
            A0(j,j) += 1;
            // Fill the off diagonal (negative part in \Delta U_{i,j}, the matrix is kept symmetric
            A0(i,j) -= 1;
            A0(j,i) -= 1;
        }

        // Fill in the rest (source and ground), i.e. components with $\Delta U_{i,j}$ with j > 15
        // Each node i connected to ground or source contributes to the rhs with R * U (R resistence between node i and ground/source
        // node, U is voltage at sink or source) and to its own diagonal with R
        rhs = Eigen::Matrix<double, nnodes, 1>::Zero();
        for(auto it: S) {
            // Shift index down by 1 and get voltage in V
            int i = std::get<0>(it) - 1;
            auto V = std::get<2>(it);
            // Add voltage in rhs (resistance assumed to be R)
            rhs(i) += V;
            // Add resistance to matrix diagonal: contribution of source current
            A0(i,i) += 1;
        }
        
        // Precompute lu factorizaion of A0
        lu = A0.lu();
    };
    
    //! \brief Compute the impedance given the resistence variable R_x
    //! Use SMW formula for low rank perturbations to avoid expensive inversion
    //! if factorization of the base system matrix is already known
    //! \param Rx resistence R_x > 0 of the varistor between node 14 and 15
    //! \return impedance = U / V of the system A(R_x)
    double operator()(double Rx) {
        // Store the scaled factor for convenience
        double fac = R/Rx;
        
        // There are many ways to create the same matrix u*v
        // Create u, the 15x2 matrix in (A+UV)
        Eigen::Matrix<double, nnodes, 2> u = Eigen::Matrix<double, nnodes, 2>::Zero();
        u(13,0) = 1;
        u(14,0) = -1;
        u(13,1) = -1;
        u(14,1) = 1;
        // v will be u tranpsosed * (1-fac)
        
        //// Use SMW formula to compute (A + UV)^{-1} rhs
        // Formula is A^{-1} rhs - A^{-1} * u * (I + V*A^{-1}*U)^{-1} V * A^{-1}
        
        // Start by precomputing A^{-1} rhs (needed twice), column vector of length 15
        auto Ainvrhs = lu.solve(rhs);
        // The precompute A^{-1} U (15x2 matrix), needed twice
        auto Ainvu = lu.solve(u);
        // The compute alpha, 2x2 matrix whose inverse is cheap
        auto alpha = (Eigen::Matrix2d::Identity() + u.transpose() / -2 * (1-fac) * Ainvu).inverse();
        // Put the formula toghether, x is a 15x1 column vector containing voltages at each node (except 16,17, prescribed)
        auto x = Ainvrhs - Ainvu * alpha * u.transpose() / -2 * (1-fac) * Ainvrhs;
        
        // Compute the current I = \Delta U_{16,5} / R and then impedance = U / I
        // \Delta U_{16,5} = (U - x_5)
        return U * R / (U - x(5));
    };
private:
    Eigen::PartialPivLU<Eigen::Matrix<double, nnodes, nnodes> > lu; //< Store lu decomposition of a for efficiency
    double R, U; //< Resistance R = R1 = R2 = R3 = R4 and source voltage U
    Eigen::Matrix<double, nnodes, 1> rhs; //< Store rhs vector prescribing sink and source voltages
};

int main(void) {
    ImpedanceMap IM = ImpedanceMap(1, 1);
    std::cout << "Impedance [R = 1]: " << IM(1) << std::endl;
    
    std::cout << std::setw(30) << "Impedance [Ohm]" << std::setw(30) << "R_x [Ohm]" << std::endl;
    for(auto Rx = 1; Rx <= 1024; Rx *= 2) {
        std::cout << std::setw(30) << IM(Rx)        << std::setw(30) << " " << Rx << std::endl;
    }
    
    timer<> tmr;
    tmr.start();
    ImpedanceMap IM2 = ImpedanceMap(1, 1);
    IM2(1);
    tmr.stop();
    std::cout << "Took:             " << tmr.min().count() / 1000000. << " ms" << std::endl;
}
