/* 

CFS-PML scheme for fractional viscoacoustic simulation (Amax = pi*f_0 CFS-PML, Memory Length L = 300)

Copyright (C) 2020 China University of Geosciences, Wuhan (Yufeng Wang)
This program is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this
program. If not, see https://www.gnu.org/licenses/

*/

#include "omp.h"
#include "stdio.h"
#include "stdlib.h"
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>   

#include "cwp.h" 
#include "su.h"
#include "header.h"
#include <signal.h>

using namespace std;

#include "Myfunctions.h"
#include "Myfunctions.c"

int main(int argc, char* argv[]) 
{

	int iz,ix,it,is;

	//****************************************************************
	//***************************define parameters********************
	//****************************************************************

	float f0=30.0;
	float t0=1/f0;
	float Omega0=20*pi*f0;	
	int PML=20;
	int Lc=8;
	int PMLc = PML + Lc;
	int nz0=234;
	int nx0=663;	
	int nz=nz0+2*PMLc;
	int nx=nx0+2*PMLc;
	int nxh=nx/2+1;
	int nzh=nz/2+1;		
	int nt=8000;
	float dt=0.001;
	float dz=10.0;
	float dx=10.0;
	int s_z=PML+2*Lc; 
	int r_z=PML+2*Lc; 
	int ds=1000;
	int ns=(nx-2*(PML+5))/ds+1;

	int Memory_len=300;	//>=4

	clock_t start, end;

	FILE *fp;
	char filename[40];

	//****************Space accurate**************
	float *rc;
	rc=(float*)malloc(sizeof(float)*Lc);
	cal_coefficients(Lc,rc);
	printf("The space accurate is %f-th order.\n",rc[0]);	


	int *s_x;		//define source position
	s_x=(int *) malloc(ns*sizeof(int*));
	for(is=0;is<=ns-1;is++)
		s_x[is]=PML+Lc+250;//is*ds+PML+5;

	float *ricker, *INT_ricker;	//Ricker wavelet
	ricker=(float *) malloc(nt*sizeof(float*));
	INT_ricker=(float *) malloc(nt*sizeof(float*));
	for(it=0;it<=nt-1;it++)
		ricker[it]=(1-2.0*(float)powf(pi*f0*(dt*it-t0),2))*exp(-(float)powf(pi*f0*(dt*it-t0),2));
	INT_ricker[0]=ricker[0];
	for(it=1;it<=nt-1;it++)  
	{
		INT_ricker[it]=INT_ricker[it-1]+ricker[it];
	}


	//****************************************************************
	//*****************************velocity and Q***********************
	//****************************************************************
	
	float **initial_vel, **vel, **vel_smooth, **initial_Q, **Q, **Gamma, ***Coeff, **receive, *energydecay;
	initial_vel =alloc2float(nx0,nz0); 		//define an initial velocity model
	vel =alloc2float(nx,nz); 
	vel_smooth =alloc2float(nx,nz);
	initial_Q =alloc2float(nx0,nz0); 		//define an initial Q model
	Q =alloc2float(nx,nz);
	Gamma =alloc2float(nx,nz);
	Coeff=alloc3float(nx,nz,Memory_len);
	receive=alloc2float(nt,4);
	energydecay = (float *) malloc(nt*sizeof(float*));
	for(it=0;it<=nt-1;it++)  
	{
		energydecay[it]=0.0;
	}

	fp=fopen("./input/acc_vp.dat", "rb");
	for (ix=0;ix<nx0;ix++)
		for(iz=0;iz<nz0;iz++)
			fread(&initial_vel[iz][ix],sizeof(float),1,fp);				
	fclose(fp); 

	fp=fopen("./input/acc_Qp.dat", "rb");
	for (ix=0;ix<nx0;ix++)
		for(iz=0;iz<nz0;iz++)
			fread(&initial_Q[iz][ix],sizeof(float),1,fp);				
	fclose(fp); 		


	createvel_Q(nx, nz, PMLc, initial_vel, vel);
	createvel_Q(nx, nz, PMLc, initial_Q, Q);

	fp=fopen("./output/vel.dat", "wb");
	for (ix=0;ix<nx;ix++)
		for(iz=0;iz<nz;iz++)
			fwrite(&vel[iz][ix],sizeof(float),1,fp);				
	fclose(fp); 

	fp=fopen("./output/Q.dat", "wb");
	for (ix=0;ix<nx;ix++)
		for(iz=0;iz<nz;iz++)
			fwrite(&Q[iz][ix],sizeof(float),1,fp);				
	fclose(fp); 

	for(iz=0;iz<=nz-1;iz++)  
		for(ix=0;ix<=nx-1;ix++)
			Gamma[iz][ix]=atan(1/Q[iz][ix])/pi;

	
	for(iz=0;iz<=nz-1;iz++)  
		for(ix=0;ix<=nx-1;ix++)
		{
			Coeff[0][iz][ix]=1.0;

			for(int j=1;j<=Memory_len-1;j++)
				Coeff[j][iz][ix]=(1-(1+2-2*Gamma[iz][ix])/j)*Coeff[j-1][iz][ix];
		}
	

	//****************************************************************
	//***************************CPML parameters**********************
	//****************************************************************
	float R=1E-6;
	float Vmax=3000.0;
	float Gmax=1.0;
	float Amax=1.0*pi*f0;

	float *dampx,*dampz, *gammax, *gammaz, *alphax, *alphaz;
	float *Ax, *Bx, *Cx, *Az, *Bz, *Cz;
	dampx=alloc1float(nx); gammax=alloc1float(nx); alphax=alloc1float(nx); 
	dampz=alloc1float(nz); gammaz=alloc1float(nz); alphaz=alloc1float(nz);
	Ax=alloc1float(nx); Bx=alloc1float(nx); Cx=alloc1float(nx);
	Az=alloc1float(nz); Bz=alloc1float(nz); Cz=alloc1float(nz);
	createCPMLpar(nx, nz, dt, dx, dz, R, Vmax, Gmax, Amax, PMLc,
		dampx, dampz, gammax, gammaz, alphax, alphaz, 
		Ax, Bx, Cx, Az, Bz, Cz);	



	//****************************************************************
	//***************************wavefileds****************************
	//****************************************************************

	float ***u;
	u=alloc3float(nx,nz,Memory_len);

	float **record, **record_remove;
	record=alloc2float(nx,nt); 
	record_remove=alloc2float(nx,nt);


	float **AX0,**AX1,**BX0,**BX1,**CX0,**CX1,**DX0,**DX1,**EX0,**EX1,**FX0,**FX1,**GX0,**GX1,**HX0,**HX1;
	float **AZ0,**AZ1,**BZ0,**BZ1,**CZ0,**CZ1,**DZ0,**DZ1,**EZ0,**EZ1,**FZ0,**FZ1,**GZ0,**GZ1,**HZ0,**HZ1;
		
	AX0=alloc2float(nx,nz);AX1=alloc2float(nx,nz);BX0=alloc2float(nx,nz);BX1=alloc2float(nx,nz),CX0=alloc2float(nx,nz);CX1=alloc2float(nx,nz);
	DX0=alloc2float(nx,nz);DX1=alloc2float(nx,nz);EX0=alloc2float(nx,nz);EX1=alloc2float(nx,nz);FX0=alloc2float(nx,nz);FX1=alloc2float(nx,nz);
	GX0=alloc2float(nx,nz);GX1=alloc2float(nx,nz);HX0=alloc2float(nx,nz);HX1=alloc2float(nx,nz);
	AZ0=alloc2float(nx,nz);AZ1=alloc2float(nx,nz);BZ0=alloc2float(nx,nz);BZ1=alloc2float(nx,nz);CZ0=alloc2float(nx,nz);CZ1=alloc2float(nx,nz);
	DZ0=alloc2float(nx,nz);DZ1=alloc2float(nx,nz);EZ0=alloc2float(nx,nz);EZ1=alloc2float(nx,nz);FZ0=alloc2float(nx,nz);FZ1=alloc2float(nx,nz);
	GZ0=alloc2float(nx,nz);GZ1=alloc2float(nx,nz);HZ0=alloc2float(nx,nz);HZ1=alloc2float(nx,nz);


	start=clock();

	for(is=0;is<=ns-1;is++)
	{

		//****************************************************************
		//*********attenuation forward propagation with vel ********************
		//****************************************************************

		printf("%d shots attenuation forward propagation with vel begin: \n",is+1);

		initialization(Memory_len, nx, nz, u,
			AX0,AX1,BX0,BX1,CX0,CX1,DX0,DX1,
			EX0,EX1,FX0,FX1,GX0,GX1,HX0,HX1,
			AZ0,AZ1,BZ0,BZ1,CZ0,CZ1,DZ0,DZ1,
			EZ0,EZ1,FZ0,FZ1,GZ0,GZ1,HZ0,HZ1);
	

		for(it=0;it<=nt-1;it++)  
		{
			int forwardorbackward=1;

			visco_CPML_FMP_Modeling(forwardorbackward,
				it, Lc, rc, Memory_len, nx, nz, dx, dz, dt, 
				Coeff, u, vel, Gamma, Omega0,
				dampx, dampz, gammax, gammaz, alphax, alphaz, 
				Ax, Bx, Cx, Az, Bz, Cz,
				AX0,AX1,BX0,BX1,CX0,CX1,DX0,DX1,
				EX0,EX1,FX0,FX1,GX0,GX1,HX0,HX1,
				AZ0,AZ1,BZ0,BZ1,CZ0,CZ1,DZ0,DZ1,
				EZ0,EZ1,FZ0,FZ1,GZ0,GZ1,HZ0,HZ1);

		#pragma omp parallel for private(iz,ix) 
			for(iz=0;iz<nz;iz++)
			{
				for(ix=0;ix<nx;ix++)
				{
					if(iz==s_z&&ix==s_x[is])
						u[Memory_len-1][iz][ix]+=dt*dt*ricker[it];
				}
			}			

		#pragma omp parallel for private(ix)
			for(ix=0;ix<=nx-1;ix++)
				record[it][ix]=u[Memory_len-1][r_z][ix];

			if(it>=2 && it<Memory_len)
			{

				receive[0][it]=u[it][80][50];
				receive[1][it]=u[it][160][200];
				receive[2][it]=u[it][200][450];
				receive[3][it]=u[it][30][640];

				//extend iz+370 ix+250

			}
			else if(it>=Memory_len)
			{
				receive[0][it]=u[Memory_len-1][80][50];
				receive[1][it]=u[Memory_len-1][160][200];
				receive[2][it]=u[Memory_len-1][200][450];
				receive[3][it]=u[Memory_len-1][30][640];
			}

			if(it%100==0)
			{
				sprintf(filename,"./output/visco_1_1_snapshots%d_%d.dat",is+1,it);
				fp=fopen(filename, "wb");

				for (ix=Lc;ix<nx-Lc;ix++)
					for(iz=Lc;iz<nz-Lc;iz++)
						fwrite(&u[Memory_len-1][iz][ix],sizeof(float),1,fp);										
				fclose(fp); 
			}

			if(it>=0 && it<Memory_len)
			{
				for (ix=Lc;ix<nx-Lc;ix++)
						for(iz=Lc;iz<nz-Lc;iz++)
							energydecay[it] +=  u[it][iz][ix]*u[it][iz][ix];
			}
			else if(it>=Memory_len)
			{
				for (ix=Lc;ix<nx-Lc;ix++)
						for(iz=Lc;iz<nz-Lc;iz++)
							energydecay[it] += u[Memory_len-1][iz][ix]*u[Memory_len-1][iz][ix];
			}		
										

			update(is, it, Memory_len, nx, nz, u,
				AX0,AX1,BX0,BX1,CX0,CX1,DX0,DX1,
				EX0,EX1,FX0,FX1,GX0,GX1,HX0,HX1,
				AZ0,AZ1,BZ0,BZ1,CZ0,CZ1,DZ0,DZ1,
				EZ0,EZ1,FZ0,FZ1,GZ0,GZ1,HZ0,HZ1);

			if(it%100==0)
				printf("vel time = %d \n", it);
		}
	}		

	end=clock();
	printf("The running time was: %fs. \n",(float)(end-start)/CLOCKS_PER_SEC);


	for(int r=0;r<4;r++)
	{
		sprintf(filename,"./output/visco_1_1_receive%d.dat",r);
		fp=fopen(filename, "wb");
		for(it=0;it<nt;it++)
			fwrite(&receive[r][it],sizeof(float),1,fp);     
		fclose(fp);  
	}

	sprintf(filename,"./output/visco_1_1_profile.dat");
	fp=fopen(filename, "wb");
	for(ix=PMLc;ix<=nx-1-PMLc;ix++)
		for(it=0;it<nt;it++)
			fwrite(&record[it][ix],sizeof(float),1,fp);     
	fclose(fp);  	

	sprintf(filename,"./output/energydecay.dat");
	fp=fopen(filename, "wb");
	for(it=0;it<nt;it++)
		fwrite(&energydecay[it],sizeof(float),1,fp);     
	fclose(fp);  


	//*******************************************
	//******************free*********************
	//*******************************************
	free(rc); free(ricker); free(INT_ricker); free(s_x);

	free2float(initial_vel); free2float(vel); free2float(vel_smooth);
	free2float(initial_Q); free2float(Q); free2float(Gamma); free3float(Coeff); free2float(receive); 
	free(energydecay);
	free3float(u); free2float(record); free2float(record_remove);

	free(Ax);free(Bx);free(Cx);free(Az);free(Bz);free(Cz);
	free(dampx);free(gammax);free(alphax);
	free(dampz);free(gammaz);free(alphaz);
	
	free2float(AX0);free2float(AX1);free2float(BX0);free2float(BX1);free2float(CX0);free2float(CX1);
	free2float(DX0);free2float(DX1);free2float(EX0);free2float(EX1);free2float(FX0);free2float(FX1);
	free2float(GX0);free2float(GX1);free2float(HX0);free2float(HX1);
	
	free2float(AZ0);free2float(AZ1);free2float(BZ0);free2float(BZ1);free2float(CZ0);free2float(CZ1);
	free2float(DZ0);free2float(DZ1);free2float(EZ0);free2float(EZ1);free2float(FZ0);free2float(FZ1);
	free2float(GZ0);free2float(GZ1);free2float(HZ0);free2float(HZ1);
	

	return 0;
}
