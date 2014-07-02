/*

gcc l.c -Wall -lm 
time ./a.out



numerical approximation


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// THE REAL SLICE OF THE MANDELBROT SET 
long double CxMin = -1.42L; //-1.4011551890L; // The Feigenbaum Point = the limit of the period doubling cascade of bifurcations 
long double CxMax = 0.26L;
long double dC ; //= CxMax - CxMin;  

// long double Cy = 0.0L;



// non linear scale of x axis 
long double lMax =  0.0L;  
long double lMin = -11.0L;
long double dl ; // = lMax - lMin;
unsigned int iPixelsNumber = 100; // number of points to check 
long double  step_l; //  = dl/iPixelsNumber;
//
long double eMin ; // = expl(lMin); // 0.0 < eMin < 1.0 
long double eMax ; // = expl(lMax); // 1.0
long double de ; // = eMax - eMin;
                                   

long double GiveC( long double l)
{
  return CxMin + expl(l)*dC;
}



int main()


{

long double Cx;
long double l;

// setup 
dC = CxMax - CxMin;  
dl = lMax - lMin;
step_l = dl/iPixelsNumber;
eMin = expl(lMin); // = expl(lMin); // 0.0 < eMin < 1.0 
eMax = expl(lMax); // 1.0
de = eMax - eMin;

// go along real axis from CxMin to CxMax using nonlinear scale 
l = lMin;
Cx = GiveC(l);
printf("c = %.20Lf \n",Cx);
while (l<lMax)
{ 
  // info message
  Cx = GiveC(l);
  printf("c = %.20Lf \n",Cx);
  // next  point 
  l += step_l;
}

  printf("nonlinear scale \n" );
return 0;
}
 
