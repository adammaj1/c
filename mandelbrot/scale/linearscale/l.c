/*

gcc l.c -Wall -lm 
time ./a.out



numerical approximationprintf("c = %.20Lf \n",Cx);


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>





int main()


{
// THE REAL SLICE OF THE MANDELBROT SET 
long double CxMin = -1.4011551890L; // The Feigenbaum Point = the limit of the period doubling cascade of bifurcations 
long double CxMax = 0.26L;  
long double Cx;
// long double Cy = 0.0L;
unsigned int iPixelsNumber = 100;
long double PixelWidth = (CxMax-CxMin)/iPixelsNumber;

// go along real axis from CxMin to CxMax using linear scale 
Cx = CxMin;
while (Cx<CxMax)
{ 
  // info message
  printf("c = %.20Lf \n",Cx);
  // next c point 
  Cx += PixelWidth;
}

  printf("linear scale with constant PixelWidth = %.20Lf \n",PixelWidth);
return 0;
}
 
