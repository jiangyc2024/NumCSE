#ifndef CIRCUITIMPEDANCE_HPP
#define CIRCUITIMPEDANCE_HPP

#include <vector>

#include <Eigen/Dense>
#include <Eigen/LU>

using namespace Eigen;

using resistor = std::pair<int, int>;
using resistor_topology = std::vector<resistor>;
using voltage = std::tuple<int, int, double>;
using voltage_topology = std::vector<voltage>;

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
	VectorXd nodal_voltage_;
	double R_, Rx_, xi_;
};

NodalPotentials::NodalPotentials(double R, double Rx): R_(R), Rx_(Rx) {
	assert(Rx_ != 0. && "Rx should not be equal to zero!");
	xi_ = R_ / Rx_;
	
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
	// We omit the resistor between node 14 and 15 for an easier matrix filling.
	
	// We implement voltage ground and sources by
	// specifying which node is connected to ground/source
	// trough a resistance. This will be also part of the r.h.s.
	voltage_topology S;
	S.reserve(4);
	S.push_back(voltage(6, 16, 1));
	S.push_back(voltage(7, 17, 0));
	S.push_back(voltage(11,17, 0));
	S.push_back(voltage(14,17, 0));
	
	// Automatically build A_Rx filling from topology
	// A_Rx is the matrix $A_R$ (i.e. with $R_x = R$)
	MatrixXd A_Rx = MatrixXd::Zero(15, 15);
	for(resistor & res: T) {
		// Shift indices down by 1 (indices in Figure start from 1)
		std::size_t i = res.first - 1;
		std::size_t j = res.second - 1;
		// Fill diagonal: each resistance contributes += R
		// into the diagonal of both nodes
		A_Rx(i,i) += 1;
		A_Rx(j,j) += 1;
		// Fill the off diagonal (negative part in $\Delta W_{i,j}$,
		// the matrix is kept symmetric
		A_Rx(i,j) -= 1;
		A_Rx(j,i) -= 1;
	}
	
	// Fill in the rest (source and ground), i.e.
	// components with $\Delta W_{i,j}$ with $j > 15$.
	// Each node $i$ connected to ground or source contributes
	// to to its own diagonal with $R$
	for(voltage& volt: S) {
		// Shift index down by 1 and get voltage in W2
		int i = std::get<0>(volt) - 1;
		double source_voltage = std::get<2>(volt);
		// Add resistance to matrix diagonal: contribution of source current
		// scaled by $R$.
		A_Rx(i,i) += 1;
	}
	
	// Add resistance between node 14 and 15
	A_Rx(13, 13) += xi_;
	A_Rx(13, 14) -= xi_;
	A_Rx(14, 13) -= xi_;
	A_Rx(14, 14) += xi_;
	
	// Set up r.h.s.
	VectorXd b = VectorXd::Zero(15);
	b(5) = 1.;
	
	// Solve unscaled linear system.
	nodal_voltage_ = A_Rx.lu().solve(b);
}

Eigen::VectorXd NodalPotentials::operator()(double V) const {
	// Scale the vector.
	return nodal_voltage_ * V;
}

class ImpedanceMap {
public:
	ImpedanceMap() = delete;
	ImpedanceMap(const ImpedanceMap&) = default;
	ImpedanceMap(ImpedanceMap&&) = default;
	ImpedanceMap& operator=(const ImpedanceMap&) = default;
	ImpedanceMap& operator=(ImpedanceMap&&) = default;
	~ImpedanceMap() = default;
	
	ImpedanceMap(double R);
	double operator()(double Rx) const;
private:
	PartialPivLU<MatrixXd> lu_; //< Store LU decomp. of matrix $A$.
	double R_, V_; //< Resistance $R$ and source voltage $W$.
	VectorXd b_; //< R.h.s vector prescribing sink/source voltages.
	std::size_t nnodes_;
};

ImpedanceMap::ImpedanceMap(double R): R_(R), nnodes_(15){
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
	S.push_back(voltage(6, 16, V_));
	S.push_back(voltage(7, 17, 0));
	S.push_back(voltage(11,17, 0));
	S.push_back(voltage(14,17, 0));
	
	// Automatically build A0 filling from topology
	// A0 is the matrix $A_R$ (i.e. with $R_x = R$)
	MatrixXd A0 = MatrixXd::Zero(nnodes_, nnodes_);
	for(resistor & res: T) {
		// Shift indices down by 1 (indices in Figure start from 1)
		std::size_t i = res.first - 1;
		std::size_t j = res.second - 1;
		// Fill diagonal: each resistance contributes += R
		// into the diagonal of both nodes
		A0(i,i) += 1;
		A0(j,j) += 1;
		// Fill the off diagonal (negative part in $\Delta W_{i,j}$,
		// the matrix is kept symmetric
		A0(i,j) -= 1;
		A0(j,i) -= 1;
	}
	
	// Fill in the rest (source and ground), i.e.
	// components with $\Delta W_{i,j}$ with $j > 15$.
	// Each node $i$ connected to ground or source contributes
	// to the r.h.s $b$ with $R \cdot V$
	// ($R$ is the  resistence between node $i$ and ground/source
	// node, $V$ is voltage at sink or source)
	// and to its own diagonal with $R$
	b_ = MatrixXd::Zero(nnodes_, 1);
	for(voltage & volt: S) {
		// Shift index down by 1 and get voltage in W2
		int i = std::get<0>(volt) - 1;
		double source_voltage = std::get<2>(volt);
		// Add voltage to r.h.s. (resistance assumed to be $R$)
		b_(i) += source_voltage;
		// Add resistance to matrix diagonal: contribution of source current
		// scaled by $R$.
		A0(i,i) += 1;
	}
	
	// Precompute the 'lu' factorization of A0
	lu_ = A0.lu();
}

double ImpedanceMap::operator()(double Rx) const {
	// Store the scaled factor for convenience
	double xi = R_ / Rx;
	
	// Create $u$: the vector in $A+u*v^\top$
	VectorXd u = VectorXd::Zero(nnodes_);
	u(13) = -std::sqrt(xi);
	u(14) = std::sqrt(xi);
	// Here: v = u
	
	// Use SMW formula to compute $(A + u \cdot u^\top)^{-1} \cdot b$.
	// Formula:
	// $A^{-1} * b - A^{-1}*u**(I+V*A^{-1}*U)^{-1}V*A^{-1} b$
	
	// Start by precomputing $A^{-1} b$, needed twice
	VectorXd Ainvrhs = lu_.solve(b_);
	// Then, precompute $A^{-1} u$, needed twice
	VectorXd Ainvu = lu_.solve(u);
	// Then, compute alpha. Alpha is just a number.
	double alpha = 1. + u.dot(Ainvu);
	// Put the formula toghether, x is a column vector containing voltages
	// at each node (except 16 and 17, which are prescribed)
	VectorXd x = Ainvrhs - Ainvu * u.dot(Ainvrhs) / alpha;
	
	// Compute the current $I = \Delta W_{16,5} / R$
	// and then impedance $= V / I$.
	// Here $\Delta W_{16,5} = W_{16} - x_5 = V - x_5$.
	return V_ * R_ / (V_ - x(5));
}

#endif
