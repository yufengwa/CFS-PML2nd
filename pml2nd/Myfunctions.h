#define pi 3.1416 
#define EPS 1E-4


//**************initialization and update************
void initialization(int Memory_len, int nx, int nz, float ***u,
		float **AX0, float **AX1, float **BX0, float **BX1, float **CX0, float **CX1, float **DX0, float **DX1, 
		float **EX0, float **EX1, float **FX0, float **FX1, float **GX0, float **GX1, float **HX0, float **HX1,
		float **AZ0, float **AZ1, float **BZ0, float **BZ1, float **CZ0, float **CZ1, float **DZ0, float **DZ1, 
		float **EZ0, float **EZ1, float **FZ0, float **FZ1, float **GZ0, float **GZ1, float **HZ0, float **HZ1);
void update(int is, int it, int Memory_len, int nx, int nz, float ***u,
		float **AX0, float **AX1, float **BX0, float **BX1, float **CX0, float **CX1, float **DX0, float **DX1, 
		float **EX0, float **EX1, float **FX0, float **FX1, float **GX0, float **GX1, float **HX0, float **HX1,
		float **AZ0, float **AZ1, float **BZ0, float **BZ1, float **CZ0, float **CZ1, float **DZ0, float **DZ1, 
		float **EZ0, float **EZ1, float **FZ0, float **FZ1, float **GZ0, float **GZ1, float **HZ0, float **HZ1);

void createvel_Q(int nx, int nz, int PML, float **a, float **b);

//****************Modeling*********************
void cal_coefficients(int Lc,float *rc);
void createCPMLpar(int nx, int nz, float dt, float dx, float dz, float R, float Vmax, float Gmax, float Amax, int PML,
		float *dampx, float *dampz, float *gammax, float *gammaz, float *alphax, float *alphaz, float *Ax, float *Bx, float *Cx, float *Az, float *Bz, float *Cz);
void removedirectwave(int nx, int nt, int is, int it, float dx, float dz, float dt, int *s_x, int s_z, int r_z, float t0, float vel, float **a, float **b);

void visco_CPML_FMP_Modeling(int forwardorbackward,
		int it, int Lc, float *rc, int Memory_len, int nx, int nz, float dx, float dz, float dt, 
		float ***Coeff, float ***u, float **vel, float **Gamma, float Omega0,
		float *dampx, float *dampz, float *gammax, float *gammaz, float *alphax, float *alphaz, 
		float *Ax, float *Bx, float *Cx, float *Az, float *Bz, float *Cz,
		float **AX0, float **AX1, float **BX0, float **BX1, float **CX0, float **CX1, float **DX0, float **DX1, 
		float **EX0, float **EX1, float **FX0, float **FX1, float **GX0, float **GX1, float **HX0, float **HX1,
		float **AZ0, float **AZ1, float **BZ0, float **BZ1, float **CZ0, float **CZ1, float **DZ0, float **DZ1, 
		float **EZ0, float **EZ1, float **FZ0, float **FZ1, float **GZ0, float **GZ1, float **HZ0, float **HZ1);
