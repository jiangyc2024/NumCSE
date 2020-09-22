#ifndef CIRCUITIMPEDANCE_HPP
#define CIRCUITIMPEDANCE_HPP

#include <vector>

#include <Eigen/Dense>
#include <Eigen/LU>

using namespace Eigen;

/* SAM_LISTNG_BEGIN_6 */
// Models a resistance between node i,j with i < j
// (in principle, it can handle different resistances)
using resistor = std::pair<int, int>;
// Vector containing all resistances in
// the circuit (excluded the variable one, which is modelled afterwards)
using resistor_topology = std::vector<resistor>;
// Models a ground or source of voltage, indexes store the node
// connected to the ground/source through a resistance
using voltage = std::tuple<int, int, double>;
// Vector containing a voltage object for each node connected to sink or source.
using voltage_topology = std::vector<voltage>;
/* SAM_LISTING_END_6 */

// \brief Struct implementing the topology of the circuit (cf. Figure)
/* SAM_LISTING_BEGIN_7 */
struct Topology {
	/* \brief Initializes the topologies.
	 * \param V Source voltage $V$ at node $16$ (in Volt),
	 * 			ground is set to $0V$ at node $17$
	 */
	Topology(double V) {
		// TODO: (3-5.d) (optional) Fill in the vectors to model the topology.
		// START
		
		// END
	}
	
	resistor_topology T;
	voltage_topology S;
};
/* SAM_LISTING_END_7 */

// \brief Class implementing the nodal potentials with fixed Rx
/* SAM_LISTING_BEGIN_0 */
class NodalPotentials {
public:
	NodalPotentials() = delete;
	NodalPotentials(const NodalPotentials&) = default;
	NodalPotentials(NodalPotentials&&) = default;
	NodalPotentials& operator=(const NodalPotentials&) = default;
	NodalPotentials& operator=(NodalPotentials&&) = default;
	~NodalPotentials() = default;
	
	NodalPotentials(double R, double Rx);
	[[nodiscard]] VectorXd operator()(double V) const;
private:
	// TODO: (3-5.d) Specify which variables should be stored.
	// START
	
	// END
};
/* SAM_LISTING_END_0 */

/* \brief Does necessary precomputations for the class
 * \param R Resistance (in Ohm) value of $R$
 * \param Rx Resistance (in Ohm) value of $Rx$
 */
/* SAM_LISTING_BEGIN_1 */
NodalPotentials::NodalPotentials(double R, double Rx)
	// TODO: (3-5.d) initializer list.
	// START
	
	// END
	{
	// TODO: (3-5.d) Do necessary and expensive precomputations here.
	// START
	
	// END
}
/* SAM_LISTING_END_1 */

/* \param V Source voltage $V$ at node $16$ (in Volt)
 * \return Vector of potentials at nodes
 */
/* SAM_LISTING_BEGIN_2 */
VectorXd NodalPotentials::operator()(double V) const {
	// TODO: (3-5.d) Return the values of the potential at the nodes.
	// START
	
	// dummy return
	return VectorXd::Zero(15);
	
	// END
}
/* SAM_LISTING_END_2 */

/* \brief Computes impedance of the entire circuit (between node 16 and 17) exploiting
*  the SMW formula for the inversion of low rank perturbations of a matrix A0,
*  whose factorization in known in advance
*/
/* SAM_LISTING_BEGIN_3 */
class ImpedanceMap {
public:
	ImpedanceMap() = delete;
	ImpedanceMap(const ImpedanceMap&) = default;
	ImpedanceMap(ImpedanceMap&&) = default;
	ImpedanceMap& operator=(const ImpedanceMap&) = default;
	ImpedanceMap& operator=(ImpedanceMap&&) = default;
	~ImpedanceMap() = default;
	
	ImpedanceMap(double R, double V);
	double operator()(double Rx) const;
private:
	// TODO: (3-5.g) Specify which variables should be stored.
	// START
	
	// END
	std::size_t nnodes_;
};
/* SAM_LISTING_END_3 */

/* \brief Build system matrix and r.h.s. and performs a LU decomposition
 * Does necessary and expensive precomputations
 * \param R Resistance (in Ohm) value of $R$
 * \param V Source voltage $V$ at node $16$ (in Volt),
 *          ground is set to $0V$ at node $17$
 */
/* SAM_LISTING_BEGIN_4 */
ImpedanceMap::ImpedanceMap(double R, double V):
	// TODO: (3-5.g) initializer list.
	// START
	
	// END
	nnodes_(15) {
	// TODO: (3-5.g) Do necessary and expensive precomputations here.
	// START
	
	// END
}
/* SAM_LISTING_END_4 */

/* \brief Compute the impedance given the resistance $R_x$.
 * Use SMW formula for low rank perturbations
 * \param Rx Resistence $R_x > 0$ between node 14 and 15
 * \return Impedance $V / I$ of the system $A_{R_x}$
 */
/* SAM_LISTING_BEGIN_5 */
double ImpedanceMap::operator()(double Rx) const {
	// TODO: (3-5.g) Return the impedance of the circuit given a specific Rx
	// START
	
	// dummy return
	return 0.;
	
	// END
}
/* SAM_LISITNG_END_5 */

#endif
