#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/LU>

#if SOLUTION
#include "timer.h"
#endif // SOLUTION

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

/* \brief Class implementing the topology of the circuit (cf. Figure)
 * Computes impedance of the entire circuit (between node 16 and 17) exploiting
 * the SMW formula for the inversion of low rank perturbations of a matrix A0,
 * whose factorization in known in advance
 */
/* SAM_LISTING_BEGIN_1 */
class ImpedanceMap {
    std::size_t nnodes; //< Number of nodes in the circuit
public:
    /* \brief Build system matrix and r.h.s. and perform a LU decomposition
     * The LU decomposition is stored in 'lu' and can be
     * reused in the SMW formula
     * to avoid expensive matrix solves
     * for repeated usages of the operator()
     * \param R Resistance (in Ohm) value of $R$
     * \param V Source voltage $V$ at node $16$ (in Volt),
     *          ground is set to $0V$ at node $17$
     */
    ImpedanceMap(double R, double V) : R(R), V(V) {
#if SOLUTION
        // In the following, instead of specifying directly
        // the entries of A\_0 and of r.h.s.
        // we define some auxiliary structure and automatically
        // generate the right entries for
        // the structures we specify.
        // This way: we avoid silly mistakes, avoid a very long list of matrix
        // entries have a flexible way (we may easily change the topology)
        // and automatically take care of shifting the indices
        // Of course, if you want, you can manually put each
        // entry in the matrix and r.h.s.
        // The result is the same.

        // Number of nodes
        nnodes = 15;

        // We implement the topology of the resistances by pushing each
        // Resistance between node i < j in a std::vector of
        // pair $(i,j) \in \mathbb{N}^2$
        // This way we can automatically build a symmetric
        // matrix, take care of indexing
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
        // The "base" matrix $A_0$ will have no resistance between 14 and 15

        // We implement voltage ground and sources by
        // specifying which node is connected to ground/source
        // trough a resistance. This will be also part of the r.h.s.
        voltage_topology S;
        S.reserve(4);
        S.push_back(voltage(6, 16, V));
        S.push_back(voltage(7, 17, 0));
        S.push_back(voltage(11,17, 0));
        S.push_back(voltage(14,17, 0));

        // Automatically build A0 filling from topology
        // A0 is the matrix $A_R$ (i.e. with $R_x = R$)
        MatrixXd A0 = MatrixXd::Zero(nnodes, nnodes);
        for(resistor & res: T) {
            // Shift indices down by 1 (indices in Figure start from 1)
            std::size_t i = res.first - 1;
            std::size_t j = res.second - 1;
            // Fill diagonal: each resistance contributes += R
            // into the diagonal of both nodes
            A0(i,i) += 1;
            A0(j,j) += 1;
            // Fill the off diagonal (negative part in \Delta W_{i,j},
            // the matrix is kept symmetric
            A0(i,j) -= 1;
            A0(j,i) -= 1;
        }

        // Fill in the rest (source and ground), i.e.
        // components with $\Delta W_{i,j}$ with $j > 15$.
        // Each node $i$ connected to ground or source contributes
        // to the r.h.s $b$ with $R \cdot V$
        // ($R$ is the  resistence between node $i$ and ground/source
        // node, $V$ is voltage at sink or source) and to its own diagonal with $R$
        b = MatrixXd::Zero(nnodes, 1);
        for(voltage & volt: S) {
            // Shift index down by 1 and get voltage in W2
            int i = std::get<0>(volt) - 1;
            double source_voltage = std::get<2>(volt);
            // Add voltage to r.h.s. (resistance assumed to be $R$)
            b(i) += source_voltage;
            // Add resistance to matrix diagonal: contribution of source current
            // scaled by $R$.
            A0(i,i) += 1;
        }

        // Precompute the 'lu' factorization of A0
        lu = A0.lu();
#else // TEMPLATE
        // TODO: build the matrix $A_0$.
        // Compute lu factorization of $A_0$.
        // Store LU factorization in 'lu'.
        // Compute the right hand side and store it in 'b'.
#endif // TEMPLATE
    };

    /* \brief Compute the impedance given the resistance $R_x$.
     * Use SMW formula for low rank perturbations to reuse LU
     * factorization.
     * \param Rx Resistence $R_x > 0$ between node 14 and 15
     * \return Impedance $V / I$ of the system $A_{R_x}$
     */
    double operator()(double Rx) {
#if SOLUTION
        // Store the scaled factor for convenience
        double f = R/Rx;

        // Create $u$: the vector in $A+u*v^\top$
        VectorXd u = VectorXd::Zero(nnodes);
        u(13) = -std::sqrt(f);
        u(14) = std::sqrt(f);
        // Here: v = u

        // Use SMW formula to compute $(A + u \cdot u^\top)^{-1} \cdot b$.
        // Formula:
        // $A^{-1} * b - A^{-1}*u**(I+V*A^{-1}*U)^{-1}V*A^{-1} b$

        // Start by precomputing $A^{-1} b$, needed twice
        VectorXd Ainvrhs = lu.solve(b);
        // Then, precompute A^{-1} u$, needed twice
        VectorXd Ainvu = lu.solve(u);
        // Then, compute alpha. Alpha is just a number.
        double alpha = 1 + u.dot(Ainvu);
        // Put the formula toghether, x is a column vector containing voltages
        // at each node (except 16 and 17, which are prescribed)
        VectorXd x = Ainvrhs - Ainvu * u.dot(Ainvrhs) / alpha;

        // Compute the current $I = \Delta W_{16,5} / R$
        // and then impedance $= V / I$.
        // Here $\Delta W_{16,5} = W_{16} - x_5 = V - x_5$.
        return V * R / (V - x(5));
#else // TEMPLATE
        // TODO: use SMW formula to compute the solution of $A_{R_x} x = b$
        // Compute and return the impedance of the system.
#endif // TEMPLATE
    }
private:
    PartialPivLU<MatrixXd> lu; //< Store LU decomp. of matrix $A$.
    double R, V; //< Resistance $R$ and source voltage $W$.
    VectorXd b; //< R.h.s vector prescribing sink/source voltages.
};
/* SAM_LISTING_END_1 */

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

#if SOLUTION
    std::cout << "--> Runtime measurments" << std::endl;
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
              << std::setprecision(3)
              << std::scientific
              << tmr_build.duration() << " s" << std::endl;
    std::cout << "Compute time: "
              << std::setprecision(3)
              << std::scientific
              << tmr_compute.duration() << " s" << std::endl;
#endif // SOLUTION

}
