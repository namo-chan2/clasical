		
		/*
 dy[i]dt=sum(j=0,j=nstate)coef[i][j]yin[j]


 
 */ 
#define NRANSI
//#include "nrutil.h"
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>
#define NR_END 1
#define FREE_ARG char*
double complex *vector(long nl,long nh);
void free_vector(double complex *v,long nl, long nh);
void rk4(double complex y[], double complex dydx[], int nstate, double x, double h, double complex yout[],double complex coef[nstate][nstate],
	void (*derivs)(double t, double complex yin[], double complex dydt[],int nstate,double complex coef[nstate][nstate]))
{
	int i;
	double xh,hh,h6;
	double complex *dym,*dyt,*yt;
	
	dym=vector(1,nstate);
	dyt=vector(1,nstate);
	yt=vector(1,nstate);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=0;i<=nstate-1;i++){ 
		yt[i]=y[i]+hh*dydx[i];
		(*derivs)(xh,yt,dyt,nstate,coef);
	}	
	for (i=0;i<=nstate-1;i++){ 	
		yt[i]=y[i]+hh*dyt[i];
		(*derivs)(xh,yt,dym,nstate,coef);
	}	
	for (i=0;i<=nstate-1;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt,nstate,coef);
	for (i=0;i<=nstate-1;i++){
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
//		printf("%1.15e\n",creal(yout[i]-y[i]));
	}	

	free_vector(yt,1,nstate);
	free_vector(dyt,1,nstate);
	free_vector(dym,1,nstate);
}
#undef NRANSI
/*

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>
#define NR_END 1
void rk4(double y[], double dydx[], int n, double x, double h, double yout[],void (*derivs)(double, double [], double [])){

	int i;
	double xh,hh,h6,*dym,*doyt,*yt;  
}	

*/

void nrerror(char error_text[]){
	fprintf(stderr,"Numericalrecipes run-time erroe..\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system\n");
	exit(1);
}	

double complex *vector(long nl,long nh){
	double complex *v;
	v=(double complex *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double complex)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}	

void free_vector(double complex *v,long nl, long nh){
	free((FREE_ARG) (v+nl-NR_END));
}	

