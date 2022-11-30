#include "seg.hpp"
#include <libgen.h>
#include <span>

//TODO(hgratten): test with MKL (without Eigen is too slow)

//1st stage of segmentation of grayscale image
int main(int argc, char** argv)
{
	const std::span args(argv, argc);
	const std::string path = std::string(dirname(args[0])) + "/test.bmp"; //NOLINT(concurrency-mt-unsafe)

	auto img = imread::readBMP(path);
	Eigen::MatrixXd P = imread::greyscale(img);

	seg::seg(P);

}
