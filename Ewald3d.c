#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include "param.dat"

// Contiguous memory allocation //
double **alloc_2d_double(int rows, int cols){
	int i=0;
	double *data = (double *)malloc(rows*cols*sizeof(double));
	double **array= (double **)malloc(rows*sizeof(double *));
	for (i=0; i<rows; i++)
		array[i] = &(data[cols*i]);
	
	return array;
	
}

double PME_Real(double **r_pos, double **prop, double **f_rewald, double *L, double *Ext_L){
	int i=0, j=0, k=0, nx=0, ny=0, nz=0;
	double rx=0.0, ry=0.0, rz=0.0;
	double alphafac = (2.0*alpha)/sqrt(M_PI);
	double alpha2 = alpha*alpha;
	double qi=0, qj=0, rsq=0.0, sqrt_rsq=0.0, erfc_sqrt_rsq=0.0, dcdr=0.0, e_rewald=0.0;
	double dij[dimen], Dij[dimen];
	for (i=0;i<ParticleN;i++){
		for (j=0;j<dimen;j++){
			f_rewald[i][j] = 0.0; 
		}
	}
	
	for (k=0;k<dimen;k++){
		dij[k] = 0.0; 
		Dij[k] = 0.0;
	}
	
	for (i=0;i<ParticleN-1;i++){
		qi = prop[i][0];
		for (j=i+1;j<ParticleN;j++){
			qj = prop[j][0];
			
			for (k=0;k<dimen;k++){
				dij[k] = r_pos[i][k]-r_pos[j][k];
				//dij[k] = dij[k] - L[k]*round(dij[k]/L[k]); // Minimum Image Convension [0,L]
				dij[k] = dij[k] - 2*Ext_L[k]*round(0.5*dij[k]/Ext_L[k]); // Minimum Image Convension [-L,L]
			}
			
			rx = dij[0];
			ry = dij[1];
			rz = dij[2];
			
			for(nx = -PME_nmax; nx <= PME_nmax; nx++){
				for(ny = -PME_nmax; ny <= PME_nmax; ny++){
					for(nz = -PME_nmax; nz <= PME_nmax; nz++){
						
						rsq = rx*rx + ry*ry + rz*rz;
						sqrt_rsq = sqrt(rsq);
						
						if(nx == 0 && ny == 0 && nz == 0){
							//central cell: ignore i == j
							if(i != j){
								erfc_sqrt_rsq = erfc(alpha*sqrt_rsq)/sqrt_rsq;
								dcdr = ONE_PI_EPS0*qi*qj*(alphafac*exp(-alpha2*rsq) + erfc_sqrt_rsq)/rsq;

								for (k=0;k<dimen;k++){
									f_rewald[i][k] = f_rewald[i][k] + dcdr*dij[k];
									f_rewald[j][k] = f_rewald[j][k] - dcdr*dij[k];
								}
								e_rewald += ONE_PI_EPS0*qi*qj*erfc_sqrt_rsq;
							}
							else{
								//ignore, because of same particle
							}
						}
						
						else{
							//image cells
							Dij[0] = rx + 2*nx*Ext_L[0];
							Dij[1] = ry + 2*ny*Ext_L[1];
							Dij[2] = rz + 2*nz*Ext_L[2];
							
							rsq = Dij[0]*Dij[0] + Dij[1]*Dij[1] + Dij[2]*Dij[2];
							sqrt_rsq = sqrt(rsq);
							
							erfc_sqrt_rsq = erfc(alpha*sqrt_rsq)/sqrt_rsq;
							dcdr = ONE_PI_EPS0*qi*qj*(alphafac*exp(-alpha2*rsq) + erfc_sqrt_rsq)/rsq;
							
							for (k=0;k<dimen;k++){
								f_rewald[i][k] = f_rewald[i][k] + dcdr*Dij[k];
								f_rewald[j][k] = f_rewald[j][k] - dcdr*Dij[k];
							}
							e_rewald += ONE_PI_EPS0*qi*qj*erfc_sqrt_rsq;
						}
					}
				}
			}								
		}	
	}
	
	printf("PME REwald: %lf\n", e_rewald);
	return e_rewald;
	
}

double PME_Fourier(double **r_pos, double **prop, double **f_kewald, double *L, double *Ext_L){
	int i=0, j=0, kx=0, ky=0,kz=0;
	double qi=0.0, ksq=0.0, kdotr=0.0, tmp=0.0, tmp2=0.0, e_kewald=0.0;
	double mx=0.0, my=0.0 ,mz=0.0, a=0.0 ,b=0.0, akak=0.0;
	double GAMMA = -0.25/(alpha*alpha);
	double recip = 2*M_PI*ONE_PI_EPS0/(8*Ext_L[0]*Ext_L[1]*Ext_L[2]);
	double ak[2];
	
	ak[0] = 0.0;
	ak[1] = 0.0;
	
	for (i=0;i<ParticleN;i++){
		for (j=0;j<dimen;j++){
			f_kewald[i][j] = 0.0; 
		}
	}
	
	for(kx = -PME_kmax; kx <= PME_kmax; kx++){
		mx = M_PI*kx/Ext_L[0];
		for(ky = -PME_kmax; ky <= PME_kmax; ky++){
			my = M_PI*ky/Ext_L[1];
			for(kz = -PME_kmax; kz <= PME_kmax; kz++){
				mz = M_PI*kz/Ext_L[2];
				ksq = mx*mx + my*my + mz*mz;
				
				if(ksq != 0){
					ak[0] = 0;
					ak[1] = 0;
					for(i=0;i<ParticleN;i++){
						qi = prop[i][0];
						kdotr = mx*r_pos[i][0] + my*r_pos[i][1] + mz*r_pos[i][2];
						
						ak[0] += qi*cos(kdotr);
						ak[1] -= qi*sin(kdotr);
					}
					
					a = ak[0];
					b = ak[1];
					akak = (a*a + b*b);
					tmp = recip * exp(GAMMA*ksq)/ksq;
					
					for(i=0;i<ParticleN;i++){
						kdotr = mx*r_pos[i][0] + my*r_pos[i][1] + mz*r_pos[i][2];
						qi = prop[i][0];
						tmp2 = 2*tmp*qi*(sin(kdotr) * a + cos(kdotr) * b);

						f_kewald[i][0] += tmp2 * mx;
						f_kewald[i][1] += tmp2 * my;
						f_kewald[i][2] += tmp2 * mz;
					}
					e_kewald += tmp * akak;
				}
			}
		}
	}
	
	printf("PME KEwald: %lf\n", e_kewald);
	return e_kewald;
}

double PME_Self(double **r_pos, double **prop){
	int i=0;
	double e_selfewald = 0.0;
	for(i=0;i<ParticleN;i++){
		e_selfewald += prop[i][0]*prop[i][0];
	}
	e_selfewald *= -ONE_PI_EPS0*alpha/sqrt(M_PI);
	
	printf("PME Self_Ewald: %lf\n", e_selfewald);
	return e_selfewald;
} 



int main(){
	
	int i=0, j=0, k=0;
	double E_rewald=0.0, E_kewald=0.0, E_selfewald=0.0;
	double **r_pos, **prop, **f_rewald, **f_kewald;
	
	// Constucting Box boundaries //
	double L[dimen];
	double Ext_L[dimen];

	L[0]=Lx;
	L[1]=Ly;
	L[2]=Lz;
	Ext_L[0]=Ext_Lx;
	Ext_L[1]=Ext_Ly;
	Ext_L[2]=Ext_Lz;
	
	// Storing the particle's initial positions //
	r_pos = alloc_2d_double(ParticleN,dimen);
	f_rewald = alloc_2d_double(ParticleN,dimen);
	f_kewald = alloc_2d_double(ParticleN,dimen);
	
	for (i=0;i<ParticleN;i++){
		for (j=0;j<dimen;j++){
			r_pos[i][j] = 0.0;
			f_rewald[i][j] = 0.0;
			f_kewald[i][j] = 0.0;
		}
	}
	
	char filename[sizeof ParticleN];
	sprintf(filename, "%d", ParticleN);
	
	FILE *fptr;
	fptr=fopen(filename,"r");
	
	double x_pos=0.0; 
	double y_pos=0.0; 
	double z_pos=0.0;
	k=0;
	while (k<ParticleN){
		if (fscanf(fptr,"%lf %lf %lf", &x_pos, &y_pos, &z_pos) == 3){
			r_pos[k][0]=x_pos;
			r_pos[k][1]=y_pos;
			r_pos[k][2]=z_pos;
		}
		k=k+1;
	}
	
	// Print the initial positions //
	printf("================ INITIAL POSITIONS ========================= \n");
	for (i=0;i<ParticleN;i++){
		for (j=0;j<dimen;j++){
			printf("%d %lf ",i,r_pos[i][j]);
		}
		printf("\n");
	}
	
	// Read the Paricle Properties //
	prop = alloc_2d_double(ParticleN,4);
	
	FILE *fptr_prop;
	fptr_prop=fopen("particle_prop.dat","r");
	
	for (i=0;i<ParticleN;i++){
		prop[i][0] = 0.0;
		prop[i][1] = 1.0;
		prop[i][2] = default_mass;
		prop[i][3] = default_gama;
	}
	
	double CHARGE=0.0;
	double DIAMETER=0.0;
	double MASS=0.0;
	double GAMMA=0.0;
	k=0;
	while (k<ParticleN){
		if (fscanf(fptr_prop, "%lf %lf %lf %lf" , &CHARGE, &DIAMETER, &MASS, &GAMMA) == 4){
			prop[k][0]=CHARGE;
			prop[k][1]=DIAMETER;
			prop[k][2]=MASS;
			prop[k][3]=GAMMA;
		}
		k=k+1;
	}
	
	// Print the Particle Properties //
	printf("================ PARTICLE PROPERTIES ========================= \n");
	for (i=0;i<ParticleN;i++){
		printf("%d %lf %lf %lf %lf\n",i,prop[i][0],prop[i][1],prop[i][2],prop[i][3]);
	}
	printf("=======================================\n");
	
	// ====================================================================================================== //
	
	double erewald = PME_Real(r_pos,prop,f_rewald,L,Ext_L);
	double ekewald = PME_Fourier(r_pos,prop,f_kewald,L,Ext_L);
	double e_selfewald = PME_Self(r_pos,prop);
	
	printf("========= 3D Ewald Energies =========\n");
	printf("Total Ewald Energy: %lf, Madelung Constant: %lf\n", erewald+ekewald+e_selfewald, 2*(erewald+ekewald+e_selfewald)/ParticleN);
	
	printf("\n");
	for (i=0;i<ParticleN;i++){
		printf("=======================================\n");
		printf("Ewald Forces on Particle: %d\n", i);
		printf("f_Rewald: %lf %lf %lf\n",f_rewald[i][0],f_rewald[i][1],f_rewald[i][2]);
		printf("f_Kewald: %lf %lf %lf\n",f_kewald[i][0],f_kewald[i][1],f_kewald[i][2]);
		printf("f_Ewald: %lf %lf %lf\n",f_rewald[i][0]+f_kewald[i][0],f_rewald[i][1]+f_kewald[i][1],f_rewald[i][2]+f_kewald[i][2]);
		printf("=======================================\n");
	}
	
	
	return 0;
	
}


















