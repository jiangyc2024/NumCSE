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
	double xi_;
};

NodalPotentials::NodalPotentials(double R, double Rx): xi_(R / Rx) {
	
}

Eigen::VectorXd NodalPotentials::operator()(double V) const {
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
	
};

#endif
