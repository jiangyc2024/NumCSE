# FLAGS FROM: https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/ 
# Cmake file didn't work :-(
# Only works with MKL installed and all variables set
# tested with ubuntu 14.04 LTS


COMPILER = g++
#COMPILER = icc
FLAGS_SEQUENTIAL = -std=c++11  -m64 -I${MKLROOT}/include -O3 -DNDEBUG
FLAGS_PARALLEL = -std=c++11 -fopenmp -m64 -I${MKLROOT}/include -O3 -DNDEBUG
# Intel(R) MKL 11.3.2, Linux, None, GNU C/C++, Intel(R) 64, Static, LP64, Sequential)
FLAGS_LINK_SEQUENTIAL = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
# Intel(R) MKL 11.3.2, Linux, None, GNU C/C++, Intel(R) 64, Static, LP64, OpenMP threading, GNU (libomp)
FLAGS_LINK_PARALLEL =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl

all: gausstiming.eps

gausstiming.eps: mkl_sequential.dat mkl_parallel.dat nomkl_gauss.dat plotting.m
	matlab -nosplash -r "plotting; exit"

mkl_sequential.dat: mkl_sequential
	./$< > $<.dat
mkl_parallel.dat: mkl_parallel
	./$< > $<.dat
nomkl_gauss.dat: nomkl_gauss
	./$< > $<.dat
	
#nomkl_gauss.dat: nomkl_gauss
#	./nomkl_gauss > nomkl_gauss.dat	
	
mkl_sequential: main.cpp gausstiming.hpp
	$(COMPILER) $(FLAGS_SEQUENTIAL) -DEIGEN_USE_MKL_ALL $< -o $@ $(FLAGS_LINK_SEQUENTIAL)
	
mkl_parallel: main.cpp gausstiming.hpp
	$(COMPILER) $(FLAGS_PARALLEL) -DEIGEN_USE_MKL_ALL $< -o $@ $(FLAGS_LINK_PARALLEL)
	
nomkl_gauss: main.cpp gausstiming.hpp
	$(COMPILER) $(FLAGS_SEQUENTIAL) $< -o $@


clean:
	rm -f mkl nomkl_gauss
