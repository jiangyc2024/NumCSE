#include <Eigen/Dense>


void sinetransform(const Eigen::VectorXd &y, Eigen::VectorXd& c)
{
	int n = y.rows() + 1;
	std::complex<double> i(0,1);

	// prepare sinus terms
	Eigen::VectorXid x = Eigen::VectorXid::LinSpaced(n-1, 1, n);
	Eigen::VectorXd sinevals = x.unaryExpr([](std::complex<double> z){ return std::pow(std::exp(i*M_PI/n), z); }).imag();

	// transform coefficients
	Eigen::VectorXd yt(n);
	yt(0) = 0;
	yt.tail(n) = sinevals.array() * (y + y.reverse()).array(); 
	yt.tail(n) += 0.5*(y-y.reverse()); 

	sinevals = imag(exp(i*pi/n).^(1:n-1));
	yt = [0 (sinevals.*(y+y(end:-1:1)) + 0.5*(y-y(end:-1:1)))];
	c = fftreal(yt);
	s(1) = dot(sinevals,y);
	for k=2:N-1
		if (mod(k,2) == 0), s(k) = -imag(c(k/2+1));
		else, s(k) = s(k-2) + real(c((k-1)/2+1)); 
	end
end
}
