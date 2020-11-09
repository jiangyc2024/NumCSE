# Makefile for blas-example MAC/Linux/MKL
#CC=icpc
CC=g++

# FLAGS for MAC (tested w/ macOS 10.15.5):
#flags=-O3 -D_MAC_OS -framework Accelerate

# FLAGS for UBUNTU (tested w/ ubuntu 18.04.5 LTS):
#Get the BLAS package with: apt-get install libblas-dev
flags=-O3 -D_LINUX -lblas

#FLAGS for FEDORA (tested w/ Fedora 33):
#Get the BLAS package with: dnf install openblas-devel
# flags=-O3 -D_LINUX -lopenblas

# FLAGS for REDHAT (Student PC Machines):
# flags=-O3 -D_LINUX -lcblas -lblas -lg2c -L /usr/lib/atlas -latlas

# FLAGS if you use Intel MKL (Linux) (tested w/ ubuntu 18.04.5 LTS):
#flags= -O3 -D_USE_MKL -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -I /opt/intel/mkl/include -L /opt/intel/mkl/lib/intel64

# FLAGS if you use Intel MKL (MacOS) (tested w/ macOS 10.15.5):
#flags= -O3 -D_USE_MKL -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -I /opt/intel/mkl/include -L /opt/intel/mkl/lib/

all: mult_g++O3

Test:	blas1.cpp
	$(CC) -o blastest blas1.cpp $(flags) 

Sources:=BlasTimings.cpp ColumnMajorMatrix.cpp ColumnMajorMatrix_Multiplication_dot.cpp ColumnMajorMatrix_Multiplication_gemm.cpp ColumnMajorMatrix_Multiplication_gemv.cpp ColumnMajorMatrix_Multiplication_noblas.cpp

mult_g++O3: $(Sources)
	$(CC) -o mult_g++O3 $(Sources) $(flags) 