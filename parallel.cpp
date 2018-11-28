#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<array>
#include<omp.h>
using namespace std; //Avoids writng std:: everytime

#define NUM_THREADS 4
#define PI 3.1415926535897; // Value of pi
const int ns = 1000; // Number of steps

struct Y_P {
	float y; 
	float qwert; // y'
};

// User defined Functions
array<Y_P,ns+1> RK2(float init[],const array<float,ns+1> &x);
void writeintofile (string s,const array<float,ns+1> &x,const array<Y_P,ns+1> &y);

// Main 
int main()
{
	// Solving Harmonic Oscillator's EOM Using Euler
	 // Initial Conditions
	float t0 = 0, tmax = 10*PI; // Time interval for solving EOM
	
	float step = abs(tmax-t0)/ns; // step size
	array <float,ns+1> t,E_rk2;
	array <Y_P,ns+1> y[4];
	for (int i = 0;i<= ns; i++) {
		t.at(i)= i*step;		// array with steps
	}
	#pragma omp parallel num_threads(NUM_THREADS)
	{
		int ID = omp_get_thread_num();
		float init[2] = {float(ID/4.0),float(1-ID/4.0)};	
		array <Y_P,ns+1> y_rk2 =  RK2(init,t);
		#pragma omp critical
		y[ID] = y_rk2;
	}	
	writeintofile("0.dat",t,y[0]);
	 writeintofile("1.dat",t,y[1]);
	 writeintofile("2.dat",t,y[2]);
	 writeintofile("3.dat",t,y[3]);
}

array<Y_P,ns+1> RK2(float init[],const array<float,ns+1> &x){
	array <Y_P,ns+1> y;
	y.at(0).y = init[0];
	y.at(0).qwert = init[1];
	float h = x.at(1) - x.at(0);
	for (auto i = 0;i< (x.size()-1);i++){
		y[i+1].y  = y[i].y +  h*(y[i].qwert - h/2*y[i].y) ;
		y[i+1].qwert  = y[i].qwert -  h*(y[i].y + h/2*y[i].qwert)  ;
	}
	return y;
 }

void writeintofile (string s,const array<float,ns+1> &x,const array<Y_P,ns+1> &y){
	ofstream outputFile;
	outputFile.open(s);

	for (int i = 0; i < y.size(); i++)
	{
		outputFile << i << "\t" << x[i] << "\t" << y[i].y << "\t" << y[i].qwert  << endl;
	}
	
	outputFile.close();
}