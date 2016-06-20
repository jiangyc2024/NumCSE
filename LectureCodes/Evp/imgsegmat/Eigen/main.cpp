#include <Eigen/Dense>
#include <iostream>
#include <libgen.h>
//#include "imgsegmat.hpp"
#include "imread.hpp"


int main(int arc, char** argv)
{
	std::string path = std::string(dirname(argv[0])) +	"/eth3.bmp";
	auto img = readBMP(path);
	
	for (int i=0; i <3; ++i)
	{
		Eigen::MatrixXd G = getcolor(img,i);

		std::cout << G << std::endl;
		std::cout << G.rows() << " " << G.cols() << std::endl;
	}

	Eigen::MatrixXd G = greyscale(img);

	std::cout << G << std::endl;
	std::cout << G.rows() << " " << G.cols() << std::endl;
}
