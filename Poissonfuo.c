#include<stdio.h>									
#include<math.h>
#include<complex.h>	//This library is declared before fftw3.h
#include<fftw3.h>	// FFTW Library for DFT

#define PI 3.14159265358	//define value of pi

void writeintofile(char *p, double* x,double* y);//function for 
//writing into file
const int Npoints = 101; // grid size you want
int main(void)
{
	int i;	// iteration variable
	float c = 1; // Constant 4*pi*G/h*h
	fftw_complex *rho_in, *rho_f_out,*phi_f_in,*phi_out;
	//Complex dynamic variables f suggest fourier coefficents 
	fftw_plan fourier_transform,invfourier_transform;
	//plans for FFT and IFFT
	double *rho_data, *phi_data;
	//for writing into file
	//Memory Allocation
	rho_data = (double*) fftw_malloc(sizeof(double)*Npoints);
	phi_data = (double*) fftw_malloc(sizeof(double)*Npoints);
	rho_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);			//memory allocation
	rho_f_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);		//memory allocation
	phi_f_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);	    //memory allocation
	phi_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);		//memory allocation
	fourier_transform = fftw_plan_dft_1d(Npoints, rho_in, rho_f_out, FFTW_FORWARD, FFTW_ESTIMATE); 
	invfourier_transform = fftw_plan_dft_1d(Npoints, phi_f_in, phi_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	printf("\nInput for mass density:\n\n");
	for(i = 0; i < Npoints; i++)
	{
		if ( i == 50 )
			rho_in[i] = 1;
			//Point masses at given positions
		else	
			rho_in[i] = 0;
			// none otherwise
		printf("%d %11.7f %11.7f\n", i, creal(rho_in[i]), cimag(rho_in[i]));		//creal and cimag are functions of complex.h 
	}
printf("\n");
fftw_execute(fourier_transform); 								//Execution of FFT
printf("Fourier Transform of rho:\n\n");
	for(i = 0; i < Npoints; i++)
	{
		printf("%d %11.7f %11.7f\n", i, creal(rho_f_out[i]), cimag(rho_f_out[i]));

	}
	for(i = 0; i < Npoints; i++)
	{
		if ( i == 0)
			phi_f_in[i] = 0;
			// Setting DC component as 0
		else
			phi_f_in[i] = c*creal(rho_f_out[i])/(2*cos(2*i*PI/Npoints)-2) + I*(cimag(rho_f_out[i])/(2*cos(2*i*PI/Npoints)-2));
	}
	printf("\nTaking inverse fourier Transform of phi:\n\n Input: \n");
	for(i = 0; i < Npoints; i++)
	{
		printf("%d %11.7f %11.7f\n", i, creal(phi_f_in[i]), cimag(phi_f_in[i]));		//creal and cimag are functions of complex.h 
	}
	printf("\n");
	fftw_execute(invfourier_transform); 								//Execution of FFT
	printf("Output:\n\n");
	for(i = 0; i < Npoints; i++)
	{	phi_out[i] = phi_out[i]/Npoints;
		printf("%d %11.7f %11.7f\n", i, creal(phi_out[i]), cimag(phi_out[i]));
		rho_data[i] = creal(rho_in[i]);
		phi_data[i] = creal(phi_out[i])- creal(phi_out[0]) ;
	}
	// Freeing memory 
	writeintofile("fourier_data.dat",rho_data,phi_data);
	fftw_destroy_plan(fourier_transform);
	fftw_destroy_plan(invfourier_transform);			//Destroy plan
	fftw_free(rho_in);
	fftw_free(phi_f_in);			 					//Free memory
	fftw_free(rho_f_out);
	fftw_free(phi_out);			 						//Free memory
	return 0;
}
void writeintofile(char *p, double* x,double*y){
	FILE *infile;
	int i;
	infile = fopen(p,"w");
	for (i=0;i<Npoints;i++) 
	{
		fprintf(infile,"%d\t%f\t%f\n",i+1,x[i],y[i]);
	} 
	fclose(infile);
}