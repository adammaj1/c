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
long double tMax =  M_PI_2 ;  // 1/tanl(lMax) = 0.0
long double tMin = M_PI_4; // 1/logl(lMin) = 1.0
long double dt ; // = lMax - lMin;
unsigned int iPixelsNumber = 100; // number of points to check 
long double  step_t; //  = dl/iPixelsNumber;


                                   

long double GiveC( long double t)
{
  return CxMax - dC/tanl(t);
}



int main()


{

long double Cx;
long double t;

// setup 
dC = CxMax - CxMin;  
dt = tMax - tMin;
step_t = dt/iPixelsNumber;


// go along real axis from CxMin to CxMax using nonlinear scale 
t = tMin;
Cx = GiveC(t);
printf("c = %.20Lf ; t =  %.20Lf \n",Cx, t );
while (t<tMax)
{ 
  // info message
  Cx = GiveC(t);
  printf("c = %.20Lf ; t =  %.20Lf \n",Cx, t );
  // next  point 
  t += step_t;
}

  printf("nonlinear scale from tMin =  %.20Lf to tMax =  %.20Lf using step_t = %.20Lf \n", tMin , tMax, step_t );
return 0;
}
 
