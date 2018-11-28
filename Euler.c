#include<stdio.h>
#include<math.h>
int main()
{
  int i;
  float X[2500],P[2500],E[2500],Err[2500],T[2500],h;
  X[0]=0;
  P[0]=1;   
  E[0]=0.5;
  T[0]=0; 
h=2*M_PI/1000;
  FILE *fptr;
  fptr=fopen("euler.txt","w");
 for (i=0;i<=2000;i++)
     { 
      X[i+1]=X[i]+h*P[i];
      P[i+1]=P[i]-h*X[i];
      T[i+1]=T[i]+h;
      E[i+1]=0.5*(P[i+1]*P[i+1]+X[i+1]*X[i+1]); 
      Err[i+1]=(E[0]-E[i+1])/E[0];
      fprintf(fptr,"%2.4f %2.4f %2.4f %2.4f %2.4f\n",X[i],P[i],T[i],E[i],Err[i]);
      }
    }