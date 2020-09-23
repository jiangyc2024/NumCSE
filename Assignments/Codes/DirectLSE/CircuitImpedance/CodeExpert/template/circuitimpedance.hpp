#ifndef CIRCUITIMPEDANCE_HPP
#define CIRCUITIMPEDANCE_HPP

#include <Eigen/Dense>
#include <Eigen/LU>
#include <vector>

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

/* \brief Struct implementing the topology of the circuit (cf. Figure).
 * Used to build the matrix of the problem.
 * Note that the resistor between node 14 and 15 is omitted.
 * Members: T a resistor_topology, i.e. a vector of index-pairs; each resistor
 * 			has the same constant resistance; that is the reason why
 *the resistor between 14-15 is omitted S a voltage_topology, i.e. a vector of
 *tuples; (i,j,V) voltage V is applied to resistor between node i and j
 */
/* SAM_LISTING_BEGIN_7 */
struct Topology {
  /* \brief Initializes the topologies.
   * \param V Source voltage $V$ at node $16$ (in Volt),
   * 			ground is set to $0V$ at node $17$
   */
  Topology(double V) {
    // We implement the topology of the resistances by pushing each
    // Resistance between node i < j in a std::vector of
    // pair $(i,j) \in \mathbb{N}^2$
    // This way we can automatically build a symmetric
    // matrix, take care of indexing
    // and avoid mistakes during the matrix filling
    T.reserve(23);
    T.push_back(resistor(1, 2));
    T.push_back(resistor(1, 5));
    T.push_back(resistor(2, 5));
    T.push_back(resistor(2, 3));
    T.push_back(resistor(2, 14));
    T.push_back(resistor(3, 4));
    T.push_back(resistor(3, 15));
    T.push_back(resistor(4, 6));
    T.push_back(resistor(4, 15));
    T.push_back(resistor(5, 7));
    T.push_back(resistor(5, 14));
    T.push_back(resistor(6, 9));
    T.push_back(resistor(6, 15));
    T.push_back(resistor(7, 10));
    T.push_back(resistor(7, 11));
    T.push_back(resistor(8, 9));
    T.push_back(resistor(8, 12));
    T.push_back(resistor(8, 13));
    T.push_back(resistor(8, 15));
    T.push_back(resistor(9, 13));
    T.push_back(resistor(10, 11));
    T.push_back(resistor(11, 12));
    T.push_back(resistor(12, 13));
    // We omit the resistor between node 14 and 15 for an easier matrix filling.

    // We implement voltage ground and sources by
    // specifying which node is connected to ground/source
    // trough a resistance. This will be also part of the r.h.s.
    S.reserve(4);
    S.push_back(voltage(6, 16, V));
    S.push_back(voltage(7, 17, 0));
    S.push_back(voltage(11, 17, 0));
    S.push_back(voltage(14, 17, 0));
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
  VectorXd operator()(double V) const;

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
  // TODO: (3-5.d) Build the matrix and the r.h.s. vector by either using
  //		 the topology struct or by hardcoding the matrix and vector.
  // START

  // END
}
/* SAM_LISTING_END_1 */

/* \param V Source voltage $V$ at node $16$ (in Volt)
 * \return Vector of potentials at nodes
 */
/* SAM_LISTING_BEGIN_2 */
VectorXd NodalPotentials::operator()(double V) const {
  VectorXd return_vector = VectorXd::Zero(15);

  // TODO: (3-5.d) Return the values of the potential at the nodes.
  // START

  // END

  return return_vector;
}
/* SAM_LISTING_END_2 */

/* \brief Computes impedance of the entire circuit (between node 16 and 17)
 * exploiting the SMW formula for the inversion of low rank perturbations of a
 * matrix A0, whose factorization in known in advance
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
};
/* SAM_LISTING_END_3 */

/* \brief Build system matrix and r.h.s. and performs a LU decomposition
 * Does necessary and expensive precomputations
 * \param R Resistance (in Ohm) value of $R$
 * \param V Source voltage $V$ at node $16$ (in Volt),
 *          ground is set to $0V$ at node $17$
 */
/* SAM_LISTING_BEGIN_4 */
ImpedanceMap::ImpedanceMap(double R, double V)
// TODO: (3-5.g) initializer list.
// START

// END
{
  // TODO: (3-5.g) Build the matrix and the r.h.s. vector by either using
  // 		 the topology or by hardcoding the matrix and vector. Also, do
  //		 the setup phase for the solver here or find an even smarter way
  //		 for an efficient implementation
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
  double return_value = 0.;

  // TODO: (3-5.g) Return the impedance of the circuit given a specific Rx
  // START

  // END

  return return_value;
}
/* SAM_LISITNG_END_5 */

#endif
