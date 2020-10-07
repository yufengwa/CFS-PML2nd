
//********************initialization & update******************
//********************initialization & update******************
//********************initialization & update******************

void initialization(int Memory_len, int nx, int nz, float ***u,
		float **AX0, float **AX1, float **BX0, float **BX1, float **CX0, float **CX1, float **DX0, float **DX1, 
		float **EX0, float **EX1, float **FX0, float **FX1, float **GX0, float **GX1, float **HX0, float **HX1,
		float **AZ0, float **AZ1, float **BZ0, float **BZ1, float **CZ0, float **CZ1, float **DZ0, float **DZ1, 
		float **EZ0, float **EZ1, float **FZ0, float **FZ1, float **GZ0, float **GZ1, float **HZ0, float **HZ1)
{
	memset(&u[0][0][0],0,nx*nz*Memory_len*sizeof(float));

	memset(AX0[0],0,nx*nz*sizeof(float));		memset(AX1[0],0,nx*nz*sizeof(float));
	memset(BX0[0],0,nx*nz*sizeof(float));		memset(BX1[0],0,nx*nz*sizeof(float));
	memset(CX0[0],0,nx*nz*sizeof(float));		memset(CX1[0],0,nx*nz*sizeof(float));
	memset(DX0[0],0,nx*nz*sizeof(float));		memset(DX1[0],0,nx*nz*sizeof(float));
	memset(EX0[0],0,nx*nz*sizeof(float));		memset(EX1[0],0,nx*nz*sizeof(float));
	memset(FX0[0],0,nx*nz*sizeof(float));		memset(FX1[0],0,nx*nz*sizeof(float));
	memset(GX0[0],0,nx*nz*sizeof(float));		memset(GX1[0],0,nx*nz*sizeof(float));
	memset(HX0[0],0,nx*nz*sizeof(float));		memset(HX1[0],0,nx*nz*sizeof(float));

	memset(AZ0[0],0,nx*nz*sizeof(float));		memset(AZ1[0],0,nx*nz*sizeof(float));
	memset(BZ0[0],0,nx*nz*sizeof(float));		memset(BZ1[0],0,nx*nz*sizeof(float));
	memset(CZ0[0],0,nx*nz*sizeof(float));		memset(CZ1[0],0,nx*nz*sizeof(float));
	memset(DZ0[0],0,nx*nz*sizeof(float));		memset(DZ1[0],0,nx*nz*sizeof(float));
	memset(EZ0[0],0,nx*nz*sizeof(float));		memset(EZ1[0],0,nx*nz*sizeof(float));
	memset(FZ0[0],0,nx*nz*sizeof(float));		memset(FZ1[0],0,nx*nz*sizeof(float));
	memset(GZ0[0],0,nx*nz*sizeof(float));		memset(GZ1[0],0,nx*nz*sizeof(float));
	memset(HZ0[0],0,nx*nz*sizeof(float));		memset(HZ1[0],0,nx*nz*sizeof(float));
}

void update(int is, int it, int Memory_len, int nx, int nz, float ***u,
		float **AX0, float **AX1, float **BX0, float **BX1, float **CX0, float **CX1, float **DX0, float **DX1, 
		float **EX0, float **EX1, float **FX0, float **FX1, float **GX0, float **GX1, float **HX0, float **HX1,
		float **AZ0, float **AZ1, float **BZ0, float **BZ1, float **CZ0, float **CZ1, float **DZ0, float **DZ1, 
		float **EZ0, float **EZ1, float **FZ0, float **FZ1, float **GZ0, float **GZ1, float **HZ0, float **HZ1)
{
	int iz, ix;

	for(int j=0;j<Memory_len-1;j++)
	{
		memcpy(&u[j][0][0],&u[j+1][0][0],nx*nz*sizeof(float));
	}
	
	memcpy(AX0[0],AX1[0],nx*nz*sizeof(float)); memcpy(BX0[0],BX1[0],nx*nz*sizeof(float));
	memcpy(CX0[0],CX1[0],nx*nz*sizeof(float)); memcpy(DX0[0],DX1[0],nx*nz*sizeof(float));
	memcpy(EX0[0],EX1[0],nx*nz*sizeof(float)); memcpy(FX0[0],FX1[0],nx*nz*sizeof(float));
	memcpy(GX0[0],GX1[0],nx*nz*sizeof(float)); memcpy(HX0[0],HX1[0],nx*nz*sizeof(float));

	memcpy(AZ0[0],AZ1[0],nx*nz*sizeof(float)); memcpy(BZ0[0],BZ1[0],nx*nz*sizeof(float));
	memcpy(CZ0[0],CZ1[0],nx*nz*sizeof(float)); memcpy(DZ0[0],DZ1[0],nx*nz*sizeof(float));
	memcpy(EZ0[0],EZ1[0],nx*nz*sizeof(float)); memcpy(FZ0[0],FZ1[0],nx*nz*sizeof(float));
	memcpy(GZ0[0],GZ1[0],nx*nz*sizeof(float)); memcpy(HZ0[0],HZ1[0],nx*nz*sizeof(float));
}


void createvel_Q(int nx, int nz, int PML, float **a, float **b)
{
	FILE *fp1, *fp2;
	int iz, ix;

	for(ix=PML;ix<=nx-PML-1;ix++)
		for(iz=PML;iz<=nz-PML-1;iz++)
		{
			b[iz][ix]=a[iz-PML][ix-PML];
		}

	for (ix=PML;ix<=nx-PML-1;ix++)
	{
		for(iz=0;iz<=PML-1;iz++)
		{
			b[iz][ix]=b[PML][ix];
		}

		for(iz=nz-PML;iz<=nz-1;iz++)
		{
			b[iz][ix]=b[nz-PML-1][ix];
		}		
	}

	for (iz=0;iz<=nz-1;iz++)
	{
		for(ix=0;ix<=PML-1;ix++)
		{
			b[iz][ix]=b[iz][PML];
		}

		for(ix=nx-PML;ix<=nx-1;ix++)
		{
			b[iz][ix]=b[iz][nx-PML-1];
		}		
	}
}



//***********************Modeling*********************
//***********************Modeling*********************
//***********************Modeling*********************

void cal_coefficients(int Lc,float *rc)
{
	int m,i;
	float s1,s2;
	for(m=1;m<=Lc;m++)
	{
		s1=1.0;s2=1.0;
		for(i=1;i<m;i++)
		{
			s1=s1*i*i;
			s2=s2*(m*m-i*i);
		}
		for(i=m+1;i<=Lc;i++)
		{
			s1=s1*i*i;
			s2=s2*(i*i-m*m);
		}
		rc[m-1]=powf(-1.0,m+1)*s1/(s2*(2.0*m*m));
	}
}

void createCPMLpar(int nx, int nz, float dt, float dx, float dz, float R, float Vmax, float Gmax, float Amax, int PML,
		float *dampx, float *dampz, float *gammax, float *gammaz, float *alphax, float *alphaz, float *Ax, float *Bx, float *Cx, float *Az, float *Bz, float *Cz)
{
	int iz, ix;
	 for(ix=0;ix<=nx-1;ix++)
	{
		if(ix>=0&&ix<=PML-1)
		{
			dampx[ix]=(3*Vmax/float(2*PML*dx))*log(1/R)*powf(float(PML-ix)/float(PML),2);
			gammax[ix]=1+(Gmax-1)*powf(float(PML-ix)/float(PML),3);
			alphax[ix]=Amax*float(ix)/float(PML);
		}
		else if(ix>=nx-PML&&ix<=nx-1)
		{
			dampx[ix]=dampx[nx-ix-1];
			gammax[ix]=gammax[nx-ix-1];
			alphax[ix]=alphax[nx-ix-1];
		}
		else
		{
			dampx[ix]=0.0;
			gammax[ix]=1.0;
			alphax[ix]=0.0;
		}
		Ax[ix]=exp(-(alphax[ix]+dampx[ix]/gammax[ix])*dt);
		Bx[ix]=-0.5*(dampx[ix]/powf(gammax[ix],2))*dt;
		Cx[ix]=exp(-alphax[ix]*dt);
	}

	for(iz=0;iz<=nz-1;iz++)
	{
		if(iz>=0&&iz<=PML-1)
		{
			dampz[iz]=(3*Vmax/float(2*PML*dz))*log(1/R)*powf(float(PML-iz)/float(PML),2);
			gammaz[iz]=1+(Gmax-1)*powf(float(PML-iz)/float(PML),3);
			alphaz[iz]=Amax*float(iz)/float(PML);
		}
		else if(iz>=nz-PML&&iz<=nz-1)
		{
			dampz[iz]=dampz[nz-iz-1];
			gammaz[iz]=gammaz[nz-iz-1];
			alphaz[iz]=alphaz[nz-iz-1];
		}
		else
		{
			dampz[iz]=0.0;
			gammaz[iz]=1.0;
			alphaz[iz]=0.0;
		}
		Az[iz]=exp(-(alphaz[iz]+dampz[iz]/gammaz[iz])*dt);
		Bz[iz]=-0.5*(dampz[iz]/powf(gammaz[iz],2))*dt;
		Cz[iz]=exp(-alphaz[iz]*dt);
	}
}

void removedirectwave(int nx, int nt, int is, int it, float dx, float dz, float dt, int *s_x, int s_z, int r_z, float t0, float vel, float **a, float **b)
{
	int ix;
	float removeline_up[nx], removeline_down[nx];

	for(ix=0; ix<nx; ix++)
	{
		removeline_up[ix]=(sqrt(powf((ix-s_x[is])*dx,2)+powf((r_z-s_z)*dz,2))/vel-4.0*t0)/dt;
		removeline_down[ix]=(sqrt(powf((ix-s_x[is])*dx,2)+powf((r_z-s_z)*dz,2))/vel+4.0*t0)/dt;

		if(it>removeline_up[ix]&&it<removeline_down[ix])
			b[it][ix]=0.0;
		else
			b[it][ix]=a[it][ix];
	}
}


void visco_CPML_FMP_Modeling(int forwardorbackward,
		int it, int Lc, float *rc, int Memory_len, int nx, int nz, float dx, float dz, float dt, 
		float ***Coeff, float ***u, float **vel, float **Gamma, float Omega0,
		float *dampx, float *dampz, float *gammax, float *gammaz, float *alphax, float *alphaz, 
		float *Ax, float *Bx, float *Cx, float *Az, float *Bz, float *Cz,
		float **AX0, float **AX1, float **BX0, float **BX1, float **CX0, float **CX1, float **DX0, float **DX1, 
		float **EX0, float **EX1, float **FX0, float **FX1, float **GX0, float **GX1, float **HX0, float **HX1,
		float **AZ0, float **AZ1, float **BZ0, float **BZ1, float **CZ0, float **CZ1, float **DZ0, float **DZ1, 
		float **EZ0, float **EZ1, float **FZ0, float **FZ1, float **GZ0, float **GZ1, float **HZ0, float **HZ1)
{
	int iz, ix;

	if(forwardorbackward==0)
	{
		dx=(-1)*dx;
		dz=(-1)*dz;		
	}

	float diffx_u0, diffx_u1, diffz_u0, diffz_u1;
	float sum_u;
	float Lap;
	float coeff;

#pragma omp parallel for private(iz,ix,sum_u,coeff,Lap, diffx_u0, diffx_u1, diffz_u0, diffz_u1)
	
	//for(iz=Lc;iz<nz-Lc;iz++)
	//{
		//for(ix=Lc;ix<nx-Lc;ix++)	
	for(iz=1;iz<=nz-2;iz++)
	{
		for(ix=1;ix<=nx-2;ix++)
		{  
			sum_u=0.0;
			for(int j=1;j<=Memory_len-1;j++)
			{
				if (j>it)
					break;

				sum_u+=Coeff[j][iz][ix]*u[Memory_len-1-j][iz][ix];
			}

			diffx_u0=0.0;
			diffx_u1=0.0;
			diffz_u0=0.0;
			diffz_u1=0.0;

			for (int ic=0; ic<Lc; ic++)
			{
				diffx_u0+=rc[ic]*(u[Memory_len-1-2][iz][ix+1+ic]-2*u[Memory_len-1-2][iz][ix]+u[Memory_len-1-2][iz][ix-1-ic])/(0.5*powf(dx,2));
				diffx_u1+=rc[ic]*(u[Memory_len-1-1][iz][ix+1+ic]-2*u[Memory_len-1-1][iz][ix]+u[Memory_len-1-1][iz][ix-1-ic])/(0.5*powf(dx,2));
				diffz_u0+=rc[ic]*(u[Memory_len-1-2][iz+1+ic][ix]-2*u[Memory_len-1-2][iz][ix]+u[Memory_len-1-2][iz-1-ic][ix])/(0.5*powf(dz,2));
				diffz_u1+=rc[ic]*(u[Memory_len-1-1][iz+1+ic][ix]-2*u[Memory_len-1-1][iz][ix]+u[Memory_len-1-1][iz-1-ic][ix])/(0.5*powf(dz,2));				
			}

			AX1[iz][ix]=Ax[ix]*AX0[iz][ix]+Bx[ix]*(Ax[ix]*diffx_u0+diffx_u1);			
			BX1[iz][ix]=Ax[ix]*BX0[iz][ix]+Bx[ix]*(Ax[ix]*AX0[iz][ix]+AX1[iz][ix]);			
			CX1[iz][ix]=Ax[ix]*CX0[iz][ix]+Bx[ix]*(Ax[ix]*0.5*(u[Memory_len-1-2][iz][ix+1]-u[Memory_len-1-2][iz][ix-1])/dx+0.5*(u[Memory_len-1-1][iz][ix+1]-u[Memory_len-1-1][iz][ix-1])/dx);			
			DX1[iz][ix]=Ax[ix]*DX0[iz][ix]+Bx[ix]*(Ax[ix]*CX0[iz][ix]+CX1[iz][ix]);			
			EX1[iz][ix]=Ax[ix]*EX0[iz][ix]+Bx[ix]*(Ax[ix]*DX0[iz][ix]+DX1[iz][ix]);			
			FX1[iz][ix]=(1.0/powf(gammax[ix],3))*0.5*(u[Memory_len-1-1][iz][ix+1]-u[Memory_len-1-1][iz][ix-1])/dx
					+(3.0/powf(gammax[ix],2))*CX1[iz][ix]
					+(3.0/gammax[ix])*DX1[iz][ix]+EX1[iz][ix];					
			GX1[iz][ix]=Cx[ix]*GX0[iz][ix]+0.5*dt*(Cx[ix]*FX0[iz][ix]+FX1[iz][ix]);			
			HX1[iz][ix]=Cx[ix]*HX0[iz][ix]+0.5*dt*(Cx[ix]*GX0[iz][ix]+GX1[iz][ix]);


			AZ1[iz][ix]=Az[iz]*AZ0[iz][ix]+Bz[iz]*(Az[iz]*diffz_u0+diffz_u1);			
			BZ1[iz][ix]=Az[iz]*BZ0[iz][ix]+Bz[iz]*(Az[iz]*AZ0[iz][ix]+AZ1[iz][ix]);			
			CZ1[iz][ix]=Az[iz]*CZ0[iz][ix]+Bz[iz]*(Az[iz]*0.5*(u[Memory_len-1-2][iz+1][ix]-u[Memory_len-1-2][iz-1][ix])/dz+0.5*(u[Memory_len-1-1][iz+1][ix]-u[Memory_len-1-1][iz-1][ix])/dz);			
			DZ1[iz][ix]=Az[iz]*DZ0[iz][ix]+Bz[iz]*(Az[iz]*CZ0[iz][ix]+CZ1[iz][ix]);		
			EZ1[iz][ix]=Az[iz]*EZ0[iz][ix]+Bz[iz]*(Az[iz]*DZ0[iz][ix]+DZ1[iz][ix]);
			FZ1[iz][ix]=(1.0/powf(gammaz[iz],3))*0.5*(u[Memory_len-1-1][iz+1][ix]-u[Memory_len-1-1][iz-1][ix])/dz
					+(3.0/powf(gammaz[iz],2))*CZ1[iz][ix]
					+(3.0/gammaz[iz])*DZ1[iz][ix]+EZ1[iz][ix];					
			GZ1[iz][ix]=Cz[iz]*GZ0[iz][ix]+0.5*dt*(Cz[iz]*FZ0[iz][ix]+FZ1[iz][ix]);			
			HZ1[iz][ix]=Cz[iz]*HZ0[iz][ix]+0.5*dt*(Cz[iz]*GZ0[iz][ix]+GZ1[iz][ix]);

			coeff=powf(vel[iz][ix],2)*powf(cos(Gamma[iz][ix]*pi/2),2)*powf(Omega0,-2*Gamma[iz][ix])*powf(float(dt),2-2*Gamma[iz][ix]);

			Lap=	((1.0/powf(gammax[ix],2))*diffx_u1
				+(2.0/gammax[ix])*AX1[iz][ix]
				+BX1[iz][ix]
				-0.5*(gammax[ix+1]-gammax[ix-1])/dx*FX1[iz][ix]
				-0.5*(dampx[ix+1]-dampx[ix-1])/dx*GX1[iz][ix]
				+dampx[ix]*0.5*(alphax[ix+1]-alphax[ix-1])/dx*HX1[iz][ix])

				+((1.0/powf(gammaz[iz],2))*diffz_u1
				+(2.0/gammaz[iz])*AZ1[iz][ix]
				+BZ1[iz][ix]
				-0.5*(gammaz[iz+1]-gammaz[iz-1])/dz*FZ1[iz][ix]
				-0.5*(dampz[iz+1]-dampz[iz-1])/dz*GZ1[iz][ix]
				+dampz[iz]*0.5*(alphaz[iz+1]-alphaz[iz-1])/dz*HZ1[iz][ix]);

			u[Memory_len-1][iz][ix]=coeff*Lap-sum_u;
		}
	}
}
