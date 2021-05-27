#include <stdio.h>
#define NRANSI
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <complex.h>
#include <time.h>
#include "/home/nishikata/physics.h"
double complex *vector(long nl, long ndt);
void free_vector(double complex *v, long nl, long ndt);
void derivs(double t,double complex yin[],double complex dydt[],int element,double complex coef[element][element]){
	int j,l;
	for(j=0;j<=element-1;j++){
		dydt[j]=0;
	}	
	for(j=0;j<=element-1;j++){
		for(l=0;l<=element-1;l++){
			dydt[j]+=coef[j][l]*yin[l];
		}
	}	
}
void rk4(double complex yin[], double complex dydt[], int element, double t, double dt, double complex yout[],
	double complex coef[element][element], 
	void (*derivs)(double t,double complex yin[],double complex dydt[],int element,double complex coef[element][element]));

int main(void){
//initial condition 
//[scan parameter
	int i,j,k,l,m,n;
	int nspin,mgrid,element;
	double fintime,dmpara;
	FILE *fp1;
	if((fp1=fopen("para_cla.dat","r"))==NULL){
		fprintf(stderr,"ERROR1\n");
		return EXIT_FAILURE;
	}	
	fscanf(fp1,"%*s %d",&nspin);
	fscanf(fp1,"%*s %d",&mgrid);
	fscanf(fp1,"%*s %le",&fintime);
	fscanf(fp1,"%*s %le",&dmpara);
	printf("n:%d\n",nspin);
	printf("m:%d\n",mgrid);
	printf("t:%1.15e\n",fintime);
	printf("dmpara:%1.15e\n",dmpara);
	
	double exconst[nspin][nspin];
			 
	for(j=0;j<=nspin-1;j++){
		for(k=0;k<=nspin-1;k++){
			exconst[j][k]=0;
		}	
	}
	for(j=1;j<=nspin-1;j++){
		if(nspin=1){
			break;
		} else{	 
			exconst[0][j]=1.116398;//kyori=sqrt(2)/2
			exconst[j][0]=exconst[0][j];
		}	
	}


//	exconst[1]=-0.120012;//kyori sqrt(2)
//scan parameter]
	element=3*nspin;

	printf("element=%d\n",element);
	
	double complex*yin,*dydt,*yout,coef[element][element],spin[nspin][3],effmf[nspin][3],dampterm[nspin][3],consv;
	double t,dt,mf[mgrid], mfst, norm, calclationtime,Bmf[3];
	time_t time0,time1;
	FILE *fp2;
	char filepath[256];
	
	dt=fintime/mgrid;
	printf("dt:%1.15e\n",dt);
	yin=vector(1,element);
	dydt=vector(1,element);
	yout=vector(1,element);
	time0=time(NULL);
	
	sprintf(filepath,"sp%d.dat",nspin);
	if((fp2=fopen(filepath,"w"))==NULL){
		fprintf(stderr,"ERROR2\n");
		return EXIT_FAILURE;
	}	
//initial setting	
	double mfomega=2*M_PI*3.8e+12;	
	double mRyd=1.e-3*Rydberg;
	double sigma=1./sqrt(2*M_PI)*4.e-13;  //2ps 
	double t0=7.e-13;
	for(j=0;j<=mgrid-1;j++){
		t=j*dt;
		mf[j]=10*exp(-(t-t0)*(t-t0)/(2*sigma*sigma))*cos(mfomega*t);
	}	


	for(j=0;j<=nspin-1;j++){
		for(k=0;k<=3-1;k++){ //xyz
			if(k==2){
				spin[j][k]=1./2;
			}else{
				spin[j][k]=0;
			}
			yin[3*j+k]=spin[j][k];
		}
	}
		
	consv = 0;
	for(j=0;j<=nspin-1;j++){
		for(k=0;k<=3-1;k++){
			consv += spin[j][k]*spin[j][k];
		}
	}
//[exchange interaction

//exchange interaction]


	t=0;
	fprintf(fp2,"%1.15e ",t);
	for(j=0;j<=nspin-1;j++){
		fprintf(fp2,"%1.15e %1.15e %1.15e "
		,creal(spin[j][0]),creal(spin[j][1]),creal(spin[j][2]));
	}	
	fprintf(fp2,"%1.15e %1.15e \n",1./4-consv,mf[0]);

//dinamics calculation	
	for(j=0;j<=mgrid-1;j++){
		t=dt*(j);

		Bmf[0]=mf[j];
		Bmf[1]=0;
		Bmf[2]=mfomega*plank/(2*M_PI*1.65*bohr_mag);


//[coefficient matrix
		for(k=0;k<=element-1;k++){
			for(l=0;l<=element-1;l++){
				coef[k][l]=0;
			}
		}	

		for(k=0;k<=nspin-1;k++){
			for(l=0;l<=3-1;l++){
				effmf[k][l] = 0;
				dampterm[k][l] = 0;
			}
			for(l=0;l<=3-1;l++){
				for(m=0;m<=nspin-1;m++){
					effmf[k][l] += 2*exconst[k][m]*mRyd*spin[m][l];
				}
			}	
			for(l=0;l<=3-1;l++){
				effmf[k][l] += 1.65*bohr_mag*Bmf[l];
			}

			effmf[k][0] += 2*3*1.e-6*1.6e-19*spin[k][0];

			dampterm[k][0] = dmpara*(spin[k][1]*effmf[k][2]-spin[k][2]*effmf[k][1]);
			dampterm[k][1] = dmpara*(spin[k][2]*effmf[k][0]-spin[k][0]*effmf[k][2]);
			dampterm[k][2] = dmpara*(spin[k][0]*effmf[k][1]-spin[k][1]*effmf[k][0]);


			coef[3*k+0][3*k+1]=-2*M_PI/((1+dmpara*dmpara)*plank)*(effmf[k][2]+dampterm[k][2]);
			coef[3*k+0][3*k+2]=2*M_PI/((1+dmpara*dmpara)*plank)*(effmf[k][1]+dampterm[k][1]);
			coef[3*k+1][3*k+0]=2*M_PI/((1+dmpara*dmpara)*plank)*(effmf[k][2]+dampterm[k][2]);
			coef[3*k+1][3*k+2]=-2*M_PI/((1+dmpara*dmpara)*plank)*(effmf[k][0]+dampterm[k][0]);
			coef[3*k+2][3*k+0]=-2*M_PI/((1+dmpara*dmpara)*plank)*(effmf[k][1]+dampterm[k][1]);
			coef[3*k+2][3*k+1]=2*M_PI/((1+dmpara*dmpara)*plank)*(effmf[k][0]+dampterm[k][0]);
		}
		
//coefficient matrix]
//[runge=kutta
		derivs(t,yin,dydt,element,coef);
		rk4(yin,dydt,element,t,dt,yout,coef,derivs);
	
		for(k=0;k<=nspin-1;k++){
			for(l=0;l<=3-1;l++){
				spin[k][l]=yout[3*k+l];
				yin[3*k+l]=yout[3*k+l];
			}
		}
//runge=kutta]
		if(j%50==0){
			consv = 0;
			for(k=0;k<=nspin-1;k++){
				for(l=0;l<=3-1;l++){
					consv += spin[k][l]*spin[k][l];
				}
			}
			fprintf(fp2,"%1.15e ",dt*(j+1));
			for(k=0;k<=nspin-1;k++){
				fprintf(fp2,"%1.15e %1.15e %1.15e ",creal(spin[k][0]),creal(spin[k][1]),creal(spin[k][2]));
			}
			fprintf(fp2,"%1.15e %1.15e \n",1./4-creal(consv),Bmf[0]);
		}
	}	
	
	fclose(fp2);
	printf("Bz=%1.15e\n",Bmf[2]);
	printf("coefBz=%1.15e\n",1.65*bohr_mag/2*Bmf[2]);
	printf("coefBx=%1.15e\n",1.65*bohr_mag/2*10);
	printf("J=%1.15e\n",exconst[0][1]*mRyd);

	free_vector(yout,1,element);
	free_vector(dydt,1,element);
	free_vector(yin,1,element);
	time1=time(NULL);
	calclationtime=difftime(time1,time0);
	printf("calclationtime=%lf\n",calclationtime);

	return 0;

}



#undef NRANSI

