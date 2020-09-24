#include "feval_noblas.h"
//we can already do a lot of preparations in the constructor:
//the idea is to write the function as  f(x) = (uv^TA + Buv^T)x =: Mx
//where M is a Matrix that can be calculated in the constructor.
//Constructor
feval_noblas::feval_noblas(int _n, 
						   const double * _A,
						   const double * _B,
						   const double * _u,
						   const double * _v ):
n(_n),M(_n,_n)
{
	//for convenience: use ColumMajorMatrix:
	ColumnMajorMatrix A(_A,_n,_n);
	ColumnMajorMatrix B(_B,_n,_n);
	double * temp_v1=new double[n];
	double * temp_v2=new double[n];
	for (int i=0; i<n;++i)
	{
		temp_v2[i]=0;
	}
	//Multplication for temp_v1=B*u;
    B.standardVectorMultiply(_u,temp_v1);
	//Multplication for temp_v2=A^T*v;
	for(int k=0; k<n; k++)
    {
		for(int l=0; l<n; l++)
		{
			temp_v2[k] += A(l,k)*_v[l];
		}
    }
	
	for(int k=0; k<n; k++)
    {
		for(int l=0; l<n; l++)
		{
			M(k,l) = _u[k]*temp_v2[l] + temp_v1[k]*_v[l];
		}
    }
	
	delete[] temp_v1;
	delete[] temp_v2;
	
}

//evaluation function:
void feval_noblas::eval( const double * x, double * y) const
{
	assert(x!=NULL);
	assert(y!=NULL);
	M.standardVectorMultiply(x, y);
	
}
