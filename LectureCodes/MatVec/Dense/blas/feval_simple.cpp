#include "feval_simple.h"
//Most Simple Implementation 
//ignoring the fact that A,B,u,v are constant
//and not pre-computing anything in the constructor

//Constructor: allocating all the Data (also for A,B,u,v)
feval_simple::feval_simple (int _n, 
			  const double * _A,
			  const double * _B,
			  const double * _u,
			  const double * _v ):
n(_n),A(_A,_n,_n),B(_B,_n,_n){
	u=new double[n];
	v=new double[n];
	for (int i=0;i<n;++i)
	{
		u[i]=_u[i];
		v[i]=_v[i];
	}
}



void feval_simple::eval(const double * x, double * y) const
{
	double* temp_v=new double[n];
	for(int i=0;i<n;++i)
	{
		temp_v[i]=0;
	}
	double temp_s1 = 0;
	double temp_s2 = 0;
	
	for (int j=0; j<n; j++)
	{
		y[j] = 0;
		temp_s1 += v[j]*x[j];
		for(int k=0; k<n; k++)
		{
			temp_v[j] += A(j,k)*x[k];
			y[j] += B(j,k)*u[k];
		}
	}
	
	for (int k=0; k<n; k++)
	{
		y[k] = y[k]*temp_s1;
		temp_s2 += v[k]*temp_v[k];
	}
	
	for (int k=0; k<n; k++)
	{
		y[k] += temp_s2*u[k];
	}
	
	delete[] temp_v;
}