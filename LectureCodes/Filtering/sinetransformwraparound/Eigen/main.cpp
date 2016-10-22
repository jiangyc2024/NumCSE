#include <Eigen/Dense>
#include <figure.hpp>
#include "./sinetransformwraparound.hpp"

int main()
{
	int n = 20;
	Eigen::VectorXd y;
	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n, 0,2);
	y = x.array().cos();
	Eigen::VectorXd c;
	sinetransformwraparound(y, c);

	// transform back
	Eigen::VectorXd y2 = Eigen::VectorXd::Zero(n);
	for (int j=1; j<=n; ++j)
	{
		y2 += c(j-1)*(x*M_PI*j/n).array().sin().matrix();
	}

	std::cout << c << std::endl;

	mgl::Figure fig;
	fig.plot(x, y2);
	fig.title("sinetransformwraparound");
	fig.xlabel("x");
	fig.ylabel("y");
	fig.setFontSize(5);
	fig.save("sinetransformwraparound");
}
