#include <iostream>
#include "simpleTimer.h"
#include "ColumnMajorMatrix.h"
#include "feval_simple.h"
#include "feval_noblas.h"
#include "feval_blas.h"


void resetVector(int n,double * y)
{
	for (int i=0;i<n;++i)
	{
		y[i]=0;
	}
}

void initalizeData(int n, double* x, double* A, double* B,  
                   double* u, double* v){
	for (int i=0;i<n;++i)
	for (int j=0;j<n;++j)
	{
		A[n*j+i]=drand48();
		B[n*j+i]=drand48();
	}
	for (int i=0;i<n;++i)
	{
		u[i]=drand48();
		v[i]=drand48();
		x[i]=drand48();
	}
}

int main (int argc, char * const argv[]) {
	
	srand(42);
	double T0,T1,T2;
	simpleTimer watch;
	int maxPower=12;
	int rep=1;
	if (argc>1) maxPower=atoi(argv[1]);
	if (argc>2) rep=atoi(argv[2]);
	
	double err1(0),err2(0);
	double *A,*B,*u,*v,*x,*y0,*y1,*y2;
	
	printf("Timings for %i repetititions: \n", rep);
	printf("N-t(simple)-t(noblas)-t(blas)-errnoblas-errblas \n");
	
	int n=2;
	//loop for different Problem Sizes:
	for (int cpow=2;cpow<=maxPower;++cpow)
	{
		T0=T1=T2=1e6;
		n*=2;
		
		A=new double[n*n];
		B=new double[n*n];
		u = new double[n];
		v = new double[n];
		x = new double[n];
		y0 = new double[n];
		y1 = new double[n];
		y2 = new double[n];
		initalizeData(n,x,A,B,u,v);
		
		//Calling Constructor for Class feval;
		feval_simple myFeval_simple(n,A,B,u,v);
		feval_noblas myFeval_noblas(n,A,B,u,v);
		feval_blas myFeval_blas(n,A,B,u,v); 
		
		//time the Calls: (make some repetititons)
		//Calling Simple Implementation:
		for (int i=0;i<rep;++i)
		{
			resetVector(n,y0);
			watch.start();
			myFeval_simple.eval(x,y0);
			T0=min(watch.getTime(),T0);
			watch.reset();
		}
		//Calling No-Blas Implementation:
		for (int i=0;i<rep;++i)
		{
			resetVector(n,y1);
			watch.start();
			myFeval_noblas.eval(x,y1);
			T1=min(watch.getTime(),T1);
			watch.reset();
		}
		//Calling Blas Implementation:
		for (int i=0;i<rep;++i)
		{
			resetVector(n,y2);
			watch.start();
			myFeval_blas.eval(x, y2);
			T2=min(watch.getTime(),T2);
			watch.reset();
		}
		err1=0;
		err2=0;
		for (int i=0;i<n;++i)
		{
			err1+=fabs(y0[i]-y1[i]);
			err2+=fabs(y0[i]-y2[i]);
		}
		printf("%i %g %g %g %g %g \n",n,T0,T1,T2,err1,err2);
		//clean up
		delete[] A,B,u,v,x,y0,y1,y2;
	}
	
}
