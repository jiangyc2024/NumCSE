#include "prevp.hpp"
#include <libgen.h>
#include <string.h>

int main(int argc, char** argv)
{
	std::string path = std::string(dirname(argv[0])) +	"/Harvard500.mtx";
	prevp(path);
}	
