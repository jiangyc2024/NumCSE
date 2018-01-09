#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

#include <figure/figure.hpp>

// Function F from subtask b
Eigen::VectorXd F(const Eigen::Array4d &x, const Eigen::ArrayXd &t, const Eigen::ArrayXd &m) {
	double a0 = x[0], b0 = x[1], l1 = x[2], l2 = x[3];
}

// Jacobian of F, computed in subtask c
Eigen::MatrixX4d DF(const Eigen::Array4d &x, const Eigen::ArrayXd &t, const Eigen::ArrayXd &m) {
	double a0 = x[0], b0 = x[1], l1 = x[2], l2 = x[3];
}

int main() {
	std::ifstream data_file("decay.txt");
	if (!data_file.is_open()) {
		std::cerr << "Cannot open decay.txt\n";
		return 1;
	}
	std::vector<double> t_data, m_data;
	for (double t_i, m_i; data_file >> t_i >> m_i;) {
		t_data.push_back(t_i); m_data.push_back(m_i);
	}
	Eigen::Map<Eigen::VectorXd> t(t_data.data(), t_data.size()), m(m_data.data(), m_data.size());


	// Use Gauss-Newton iteration to find minimizer of F(x)
	return 0;
}
