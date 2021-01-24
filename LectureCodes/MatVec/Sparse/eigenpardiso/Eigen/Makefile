# FLAGS FROM: https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/ 
# Only works with MKL installed: https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library/choose-download.html
# MKL can be installed to /opt/intel/mkl

# Tested with ubuntu 18.04.5 LTS, eigen 3.3.7 and macOS 10.15.5, eigen 3.3.7

# MKL install location
MKLROOT = /opt/intel/mkl

# Download source here: https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip
EIGENROOT = /usr/local/include/eigen-3.3.7 

COMPILER = g++

FLAGS = -std=c++11 -DMKL_LP64 -m64 -I${MKLROOT}/include -I${EIGENROOT} -O3 -DNDEBUG 

# Ubuntu
# Intel(R) MKL 2020.0, Linux, GNU C/C++, Intel(R) 64, Static, 32-bit integer, Sequential)
#FLAGS_LINK = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
FLAGS_LINK = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -L ${MKLROOT}/lib/intel64

# MacOS
# Intel(R) MKL 2020.0, MacOS, GNU C/C++, Intel(R) 64, Static, 32-bit integer, Sequential)
#FLAGS_LINK = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -L ${MKLROOT}/lib/

all: main.cpp
	$(COMPILER) $(FLAGS) -DEIGEN_USE_MKL_ALL $< -o main $(FLAGS_LINK)

clean:
	rm -f main
