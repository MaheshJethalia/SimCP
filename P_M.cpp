#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<array>
#include <ctime>
#include <cstdlib>
#include<complex.h>	//This library is declared before fftw3.h
#include<fftw3.h>

using namespace std; //Avoids writng std:: everytime
#define PI 3.1415926535897 // Value of pi
// Number of steps
const int np = 3;	// Number of Particles
const int gns = 100; 	// Total Number of Grid Points
 
const float m = 1.0;
const float h = 0.01;

struct gh {
	float y; 
	float p; 
};

void writeintofile (string s,const array<gh,np> &y);
array<gh,np> Leapfrog(float dt, array<gh,np> &x, int Number_of_steps);
array<gh,np> generate_random_upto(int gns,float h);
array<float,gns>  rho_cic(array<gh,np> Y);
array<float,gns> del_phi(array<float,gns> &rho);

int main()
{
	float time_step = 1.0;
	int numsteps = 5; 
	array<gh,np> y,z;
	y = generate_random_upto(gns,h);
	z = Leapfrog(time_step,y,numsteps);
	return 0;
}

 array<gh,np> Leapfrog(float dt, array<gh,np> &x, int Number_of_steps){
	array<float,gns> rho ;
	array<float,gns> del_phis;
	float weighted_del_phi;	
	for(auto j = -1; j < Number_of_steps; j++){
		rho = rho_cic(x);
		del_phis = del_phi(rho);
		for (auto i = 0;i< np;i++){
			int k = floor(x[i].y/h);
			weighted_del_phi = del_phis[k]*(1 - (x[i].y/h - k)) + del_phis[k+1]*(x[i].y/h - k);
			if (j == -1){
				x[i].p = x[i].p - dt*(weighted_del_phi)/2;
			}
			else if ( j == np-1){
				x[i].p = x[i].p - dt*(weighted_del_phi)/2;
			}
			else{
				x[i].y  = x[i].y +  dt*(x[i].p)/m  ;
				if (x[i].y >= h*gns || x[i].y < 0){
					int r = floor(x[i].y/(h*gns));
					printf("\nLa x(i):%f\tL :%f",x[i].y,h*gns);
				}
				x[i].p  = x[i].p -  dt*(weighted_del_phi);
			}
		}
		writeintofile("PMfile_no" + to_string(j) + ".dat",x);
	}
	writeintofile("PMfile_no_1a.dat",x); 	
	return x;
 }

void writeintofile (string s,const array<gh,np> &y){
	ofstream outputFile;
	outputFile.open(s);

	for (int i = 0; i < y.size(); i++)
	{
		outputFile << i << "\t" << y[i].y << "\t" << y[i].p << endl;
	}
	
	outputFile.close();
}

array <gh,np> generate_random_upto(int gns,float h){
	srand(time(NULL));   
	array <gh,np> y;
	for (auto i=0;i<np;i++)
	{
		y[i].y= (float) rand()*gns*h/RAND_MAX;
		y[i].p= (float) rand()/RAND_MAX;
		printf("%f\t%f\n", y[i].y,y[i].p);
	}
	return y; 
	
}
array<float,gns>  rho_cic(array<gh,np> y){
	array <float,gns> rho;
	int i,j;
	for (i=0;i<gns;i++){
		rho[i] = 0;
	}
	for (i=0;i<np;i++){
		j = floor(y[i].y/h);
		if (j < gns && j > 0){
			rho[j] = m*(1 - (y[i].y/h - j));
			rho[j+1] = m*(y[i].y/h - j);
			} 
	}
	return rho;	
}
array<float,gns> del_phi(array<float,gns> &rho){
	int i;	
	float c = 1*h; 
	fftw_complex *rho_in, *rho_f_out,*phi_f_in,*phi_out;
	fftw_plan fourier_transform,invfourier_transform;
	array<float,gns> del_phi_data;
	//for writing into file

	//Memory Allocation			//allocating memory
	rho_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gns);			//allocating memory
	rho_f_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gns);		//allocating memory
	phi_f_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gns);			//allocating memory
	phi_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gns);		//allocating memory
	fourier_transform = fftw_plan_dft_1d(gns, rho_in, rho_f_out, FFTW_FORWARD, FFTW_ESTIMATE); 	//Here we set which kind of transformation we want to perform
	invfourier_transform = fftw_plan_dft_1d(gns, phi_f_in, phi_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	//printf("\nInput for mass density:\n\n");
	for(i = 0; i < gns; i++)
	{
		rho_in[i][0] =  rho[i];
		rho_in[i][1] =  0;
	}
	fftw_execute(fourier_transform); 								//Execution of FFT
	for(i = 0; i < gns; i++)
	{
		if ( i == 0)
			phi_f_in[i][0] = 0;
		else
			phi_f_in[i][0] = c*rho_f_out[i][0]/(2*cos(2*i*PI/gns)-2);
			phi_f_in[i][1] = c*(rho_f_out[i][1]/(2*cos(2*i*PI/gns)-2));
			phi_f_in[i][1] = phi_f_in[i][0]*sin(2*i*PI/gns)/h;
			phi_f_in[i][0] = -phi_f_in[i][1]*sin(2*i*PI/gns)/h;
	}
	fftw_execute(invfourier_transform); 								//Execution of FFT
	for(i = 0; i < gns; i++)
	{	phi_out[i][0] = phi_out[i][0]/gns;
		phi_out[i][1] = phi_out[i][1]/gns;
		del_phi_data[i] = phi_out[i][0];
	}
	fftw_destroy_plan(fourier_transform);
	fftw_destroy_plan(invfourier_transform);							//Destroy plan
	fftw_free(rho_in);
	fftw_free(phi_f_in);			 						//Free memory
	fftw_free(rho_f_out);
	fftw_free(phi_out);			 						//Free memory
	
	return del_phi_data;
}
