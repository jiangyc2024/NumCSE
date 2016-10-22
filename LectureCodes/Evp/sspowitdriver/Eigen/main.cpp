#include "sspowitdriver.hpp"

int main()
{
	int n = 6;
	Eigen::VectorXd d = Eigen::VectorXd::LinSpaced(n, 1, n);
	ssppowitdriver(d);
}
