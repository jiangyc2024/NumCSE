#include <Eigen/Dense>
#include <figure.hpp>
#include <iostream>
#include <libgen.h>
#include <mgl2/mgl.h>
#include "imread.hpp"


int main(int arc, char** argv)
{
	std::string path = std::string(dirname(argv[0])) +	"/eth.bmp";

	auto img = readBMP(path);
	Eigen::MatrixXd P = greyscale(img);

	P = P/P.maxCoeff(); // normalize date

	mglData a(P.data(), P.rows(), P.cols());
	
	mglGraph gr;
	gr.SetSize(800,800);

	gr.Title("Image plot");
	gr.Box();  
	gr.Rotate(0,90);
	gr.Tile(a, "kw");

	gr.WriteEPS("greyscale.eps");

}
