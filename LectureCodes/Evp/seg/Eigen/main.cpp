#include <libgen.h>
#include "seg.hpp"

//TODO: test with MKL (without Eigen is too slow)

//1st stage of segmentation of grayscale image
int main(int argc, char** argv)
{
	std::string path = std::string(dirname(argv[0])) +	"/test.bmp";

	auto img = readBMP(path);
	Eigen::MatrixXd P = greyscale(img);

	seg(P);

}
