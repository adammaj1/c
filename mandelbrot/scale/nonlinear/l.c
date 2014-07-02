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
long double PixelWidth = 0.3L;
long double r=4.05L; // https://en.wikipedia.org/wiki/Feigenbaum_constants
unsigned int i;

// go along real axis from CxMin to CxMax using linear scale 
Cx = CxMax;
i=0;
while (i<iPixelsNumber)
{ 
  // info message
  printf("i = %d ; c = %.20Lf ; PixelWidth =  %.20Lf \n",i,Cx, PixelWidth);
  // next c point 
  if (i % 10) PixelWidth/= r;
  Cx -= PixelWidth;
  i+=1;
}

  
return 0;
}
 
