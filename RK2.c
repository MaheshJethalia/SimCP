#include<stdio.h>
#include<math.h>
int main()
{
 int i;
 float X[2500],P[2500],E[2500],T[2500],Err[2500],h,s1,s2,p1,p2;
 X[0]=0;
 P[0]=1;
 E[0]=0.5;
 T[0]=0;
 Err[0]=0;
 h=2*M_PI/1000;
FILE *fptr;
 fptr=fopen("rk.txt","w");
 for (i=0; i<=2000; i++)
 {
  s1=P[i];
  p1=-X[i];
  s2=P[i]+h*p1;
  p2=-(X[i]+h*s2); 
  X[i+1]=X[i]+h*0.5*(s1+s2);
  P[i+1]=P[i]+h*0.5*(p1+p2);  
  T[i+1]=T[i]+h;
  E[i+1]=0.5*(P[i+1]*P[i+1]+X[i+1]*X[i+1]);
  Err[i+1]=(E[i+1]-E[0])/E[0];
  fprintf(fptr,"%2.4f %2.4f %2.4f %2.4f %2.4f\n",X[i],P[i],T[i],E[i],Err[i]);
 } 
 
}