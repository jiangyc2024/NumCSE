#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/LU>

#if SOLUTION
#include "timer.h"
#endif // SOLUTION

using namespace Eigen;

// Models a resistance between node i,j with i < j (could, in principle, handle different resistances)
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
#if SOLUTION
        // In the following, instead of specifying directly the entries of A\_0 and of r.h.s.
        // we define some auxiliary structure and automatically generate the right entries for
        // the structures we specify.
        // This way: we avoid silly mistakes, avoid a very long list of matrix
        // entries have a flexible way (we may easily change the topology)
        // and automatically take care of shifting the indices
        // Of course, if you want, you can manually put each entry in the matrix and r.h.s.
        // The result is the same.

        // Number of nodes
        nnodes = 15;

        // We implement the topology of the resistances by pushing each
        // Resistance between node i < j in a std::vector of pair $(i,j) \in \mathbb{N}^2$
        // This way we can automatically build a symmetric matrix, take care of indexing
        // and avoid mistakes during the matrix filling
        resistor_topology T;
        T.reserve(23);
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
        // Default matrix entries set to one (will be a "varistor")
        T.push_back(resistor(14,15));

        // We implement voltage ground and source topology by
        // specifying which node is connected to ground/source
        // trough a resistance. This will be also part of the r.h.s.
        voltage_topology S;
        S.reserve(4);
        S.push_back(voltage(6, 16, W));
        S.push_back(voltage(7, 17, 0));
        S.push_back(voltage(11,17, 0));
        S.push_back(voltage(14,17, 0));

        // Automatically build A0 filling from topology
        MatrixXd A0 = MatrixXd::Zero(nnodes, nnodes);
        for(resistor & it: T) {
            // Shift indices down by 1 (indices in Figure start from 1)
            std::size_t i = it.first - 1;
            std::size_t j = it.second - 1;
            // Fill diagonal: each resistance contributes += R into the diagonal of both nodes
            A0(i,i) += 1;
            A0(j,j) += 1;
            // Fill the off diagonal (negative part in \Delta W_{i,j}, the matrix is kept symmetric
            A0(i,j) -= 1;
            A0(j,i) -= 1;
        }

        // Fill in the rest (source and ground), i.e. components with $\Delta W_{i,j}$ with j > 15
        // Each node i connected to ground or source contributes
        // to the rhs with R * W (R resistence between node i and ground/source
        // node, W is voltage at sink or source) and to its own diagonal with R
        rhs = Matrix::Zero(nnodes,1);
        for(auto it: S) {
            // Shift index down by 1 and get voltage in W2
            int i = std::get<0>(it) - 1;
            auto W2 = std::get<2>(it);
            // Add voltage in rhs (resistance assumed to be R)
            rhs(i) += W2;
            // Add resistance to matrix diagonal: contribution of source current
            A0(i,i) += 1;
        }

        // Precompute the 'lu' factorization of A0
        lu = A0.lu();
#else // TEMPLATE
        // TODO:
#endif // TEMPLATE
    };

    //! \brief Compute the impedance given the resistance $R_x$
    //! Use SMW formula for low rank perturbations to avoid expensive solution of system
    //! if factorization of the base system matrix is already known
    //! \param Rx Resistence $R_x > 0$ of the varistor between node 14 and 15
    //! \return $impedance = W  * I$ of the system $A(R_x)$
    double operator()(double Rx) {
#if SOLUTION
        // Store the scaled factor for convenience
        double fac = R/Rx;

        // There are many ways to create the same matrix u*v
        // Create $U$: the 15x2 matrix in $A+UV$
        MatrixXd U = MatrixXd::Zero(nnodes, 1);
        U(13) = -1;
        U(14) = 1;
        // V will be U transposed

        // Use SMW formula to compute $(A + UU')^{-1}$ r.h.s.
        // Formula is A^{-1} rhs - A^{-1} * u * (I + V*A^{-1}*U)^{-1} V * A^{-1}

        // Start by precomputing A^{-1} rhs (needed twice), column vector of length 15
        auto Ainvrhs = lu.solve(rhs);
        // The precompute A^{-1} U (15x2 matrix), needed twice
        auto Ainvu = lu.solve(U);
        // The compute alpha, 2x2 matrix whose inverse is cheap
        auto alpha = (Matrix::Identity(2,2) + U.transpose() / -1 * (1-fac) * Ainvu).inverse();
        // Put the formula toghether, x is a 15x1 column vector containing voltages
        // at each node (except 16,17, prescribed)
        auto x = Ainvrhs - Ainvu * alpha * U.transpose() / -1 * (1-fac) * Ainvrhs;

        // Compute the current $I = \Delta W_{16,5} / R$ and then $impedance = W / I$
        // Here $\Delta W_{16,5} = (W - x_5)$
        return W * R / (W - x(5));
    };
private:
    PartialPivLU<MatrixXd> lu; //< Store LU decomposition of matrix A
    double R, W; //< Resistance R and source voltage W
    VectorXd rhs; //< Store r.h.s. vector prescribing sink and source voltages
#else // TEMPLATE
        // TODO:
#endif // TEMPLATE
};

int main(void) {
    // Create a first ImpedanceMap with resistance 1 and voltage 1
    ImpedanceMap IM = ImpedanceMap(1, 1);
    // Test Impedance at value one
    std::cout << "Impedance [R = 1]: " << IM(1) << std::endl;

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

#if SOLUTION
    Timer tmr_build, tmr_compute;
    // Measuring time it takes to build the matrix
    tmr_build.start();
    ImpedanceMap IM2 = ImpedanceMap(1, 1);
    tmr_build.stop();
    // Measuring time it takes to compute the impedance
    tmr_compute.start();
    IM2(1);
    tmr_compute.stop();
    // Output time
    std::cout << "Building time: "
              << std::scientific << std::precision(3)
              << tmr_build.duration() << " s" << std::endl;
    std::cout << "Compute time: "
              << std::scientific << std::precision(3)
              << tmr_compute.duration() << " s" << std::endl;
#endif // SOLUTION

}
