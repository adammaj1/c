/*

gcc l.c -Wall -lm 
time ./a.out



numerical approximationprintf("c = %.20Lf \n",Cx);


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// THE REAL SLICE OF THE MANDELBROT SET 
long double CxMin = -1.4011551890L; // The Feigenbaum Point = the limit of the period doubling cascade of bifurcations 
long double CxMax = 0.26L;
long double dC ; //= CxMax - CxMin;  

// long double Cy = 0.0L;



// non linear scale of x axis 
long double lMax =  M_E;  // logl(lMax) = 1.0
long double lMin = 1.0L; // logl(lMin) = 0.0
long double dl ; // = lMax - lMin;
unsigned int iPixelsNumber = 1000; // number of points to check 
long double  step_l; //  = dl/iPixelsNumber;


                                   

long double GiveC( long double l)
{
  return CxMin + logl(l)*dC;
}



int main()


{

long double Cx;
long double l;

// setup 
dC = CxMax - CxMin;  
dl = lMax - lMin;
step_l = dl/iPixelsNumber;


// go along real axis from CxMin to CxMax using nonlinear scale 
l = lMin;
Cx = GiveC(l);
printf("c = %.20Lf ; l =  %.20Lf \n",Cx, l );
while (l<lMax)
{ 
  // info message
  Cx = GiveC(l);
  printf("c = %.20Lf ; l =  %.20Lf \n",Cx, l );
  // next  point 
  l += step_l;
}

  printf("nonlinear scale from lMin =  %.20Lf to lMax =  %.20Lf using step_l = %.20Lf \n", lMin , lMax, step_l );
return 0;
}
 
