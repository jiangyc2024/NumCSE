#Makefile for blas-example MAC/Linux/MKL
#CC=icpc
CC=g++

#FLAGS for MAC:
flags=-O3 -D_MAC_OS -framework Accelerate

#FLAGS for UBUNTU:
#flags=-O3 -D_LINUX -lblas

#FLAGS for REDHAT (Student PC Machines):
# flags=-O3 -D_LINUX -lcblas -lblas -lg2c -L /usr/lib/atlas -latlas

#FLAGS if you use Intel MKL
#flags= -O3 -D_USE_MKL -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

Test:	blas1.cpp
	$(CC) $(flags) -o blastest blas1.cpp

Sources:=BlasTimings.cpp ColumnMajorMatrix.cpp ColumnMajorMatrix_Multiplication_dot.cpp ColumnMajorMatrix_Multiplication_gemm.cpp ColumnMajorMatrix_Multiplication_gemv.cpp ColumnMajorMatrix_Multiplication_noblas.cpp

mult_g++O3: $(Sources)
	$(CC) $(flags) -o mult_g++O3 $(Sources)

