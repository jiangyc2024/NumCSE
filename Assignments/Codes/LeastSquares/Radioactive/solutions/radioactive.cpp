#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

#include <figure/figure.hpp>

// Function F from subtask b
Eigen::VectorXd F(const Eigen::Array4d &x, const Eigen::ArrayXd &t, const Eigen::ArrayXd &m) {
	double a0 = x[0], b0 = x[1], l1 = x[2], l2 = x[3];
	return (-l2*t).exp()*b0 + (l1/(l2-l1))*((-l1*t).exp() - (-l2*t).exp())*a0 - m;
}

// Jacobian of F, computed in subtask c
Eigen::MatrixX4d DF(const Eigen::Array4d &x, const Eigen::ArrayXd &t, const Eigen::ArrayXd &m) {
	double a0 = x[0], b0 = x[1], l1 = x[2], l2 = x[3];
	Eigen::ArrayX4d res(t.size(), 4);
	double inv_d = 1.0 / (l2-l1);
	Eigen::ArrayXd expl1 = (-l1*t).exp();
	Eigen::ArrayXd expl2 = (-l2*t).exp();
	Eigen::ArrayXd expd = expl1-expl2;

	Eigen::VectorXd col1 = l1*inv_d*expd;
	Eigen::VectorXd col2 = expl2;
	Eigen::VectorXd col3 = l2*inv_d*inv_d*expd*a0 - l1*t*inv_d*expl1*a0;
	Eigen::VectorXd col4 = -t*expl2*b0 - l2*inv_d*inv_d*expd*a0 + l1*t*inv_d*expl2*a0;
	res << col1, col2, col3, col4;
	return res;
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
/* SAM_LISTING_BEGIN_1 */
	std::vector<double> gn_update;
	Eigen::Vector4d x(1., 1., 1., 0.1), s;
	do {
		Eigen::Vector4d s = DF(x, t, m).colPivHouseholderQr().solve(F(x, t, m));
		x -= s;
		gn_update.push_back(s.lpNorm<Eigen::Infinity>());
	} while (gn_update.back() > 1e-14);
/* SAM_LISTING_END_1 */


	mgl::Figure fig1, fig2;
	fig1.plot(t, m, "b").label("measured PhiB");
	fig1.plot(t, F(x, t, Eigen::VectorXd::Zero(t.size())), "r").label("fitted PhiB");
	Eigen::VectorXd phiA =  x[0]*(-x[2]*Eigen::ArrayXd(t)).exp();
	fig1.plot(t, phiA, "b").label("fitted PhiA");
	fig1.xlabel("t");
	fig1.legend(1, 1);
	fig1.save("solution.eps");

	fig2.setlog(false, true);
	fig2.plot(Eigen::VectorXd::LinSpaced(gn_update.size(), 1, gn_update.size()), Eigen::Map<Eigen::VectorXd>(gn_update.data(), gn_update.size()), "r");
	fig2.xlabel("iteration");
	fig2.ylabel("change");
	fig2.save("convergence.eps");
	return 0;
}
