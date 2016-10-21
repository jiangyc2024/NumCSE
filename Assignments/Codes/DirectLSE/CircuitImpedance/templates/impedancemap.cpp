#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/LU>


using namespace Eigen;

/* SAM_LISTING_BEGIN_2 */
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
/* SAM_LISTING_BEGIN_2 */

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
        // TODO: build the matrix $A_0$.
        // Compute lu factorization of $A_0$.
        // Store LU factorization in 'lu'.
        // Compute the right hand side and store it in 'b'.
    };

    /* \brief Compute the impedance given the resistance $R_x$.
     * Use SMW formula for low rank perturbations to reuse LU
     * factorization.
     * \param Rx Resistence $R_x > 0$ between node 14 and 15
     * \return Impedance $V / I$ of the system $A_{R_x}$
     */
    double operator()(double Rx) {
        // TODO: use SMW formula to compute the solution of $A_{R_x} x = b$
        // Compute and return the impedance of the system.
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
    std::cout << std::setw(30) << IM(0.1)
              << std::setw(30) << 0.1
              << std::endl;
    for(auto Rx = 1; Rx <= 1024; Rx *= 2) {
        std::cout << std::setw(30) << IM(Rx)
                  << std::setw(30) << Rx
                  << std::endl;
    }


}
