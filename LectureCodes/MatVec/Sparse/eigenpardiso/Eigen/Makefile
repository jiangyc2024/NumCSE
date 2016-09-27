# FLAGS FROM: https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/ 
# Cmake file didn't work :-(
# Only works with MKL installed and all variables set
# tested with ubuntu 14.04 LTS

COMPILER = g++
FLAGS = -std=c++11 -m64 -I${MKLROOT}/include -O3 -DNDEBUG

# Intel(R) MKL 11.3.2, Linux, None, GNU C/C++, Intel(R) 64, Static, LP64, Sequential)
FLAGS_LINK = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl

all: main.cpp
	$(COMPILER) $(FLAGS) -DEIGEN_USE_MKL_ALL $< -o main $(FLAGS_LINK)

clean:
	rm -f main
