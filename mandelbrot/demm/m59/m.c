// gcc m.c -Wall -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double pixelwidth = 0.003;
double bailout = 4.0;
unsigned int maxiter=255;

// mndlbrot.cpp from Mandel 5.9 by Wolf Jung (C) 2007-2013
// c=a+b*i z=x+y*i

int dist(double a, double b, double x, double y)
{  

   unsigned int j; 
   double xp = 1; //zp = 1
   double yp = 0;
   double nz, nzp; 

   for (j = 1; j <= maxiter; j++)
    { nz = 2*(x*xp - y*yp); 
      yp = 2*(x*yp + y*xp); 
      xp = nz; //zp = 2*z*zp;
      xp++; //zp = 2*z*zp + 1
      nz = x*x - y*y + a; 
      y = 2*x*y + b; x = nz; //z = z*z + c;
      nz = x*x + y*y; 
      nzp = xp*xp + yp*yp;
      if (nzp > 1e40 || nz > bailout) break;
    }
   if (nz < bailout) return 1; //not escaping, rare
   if (nzp < nz) return 10; //includes escaping through 0
   x = log(nz); 
   x = x*x*nz / (nzp); //4*square of dist/pixelwidth
   if (x < 0.04) return 1; 
   if (x < 0.24) return 9; // near boundary
   return 10;
} //dist


int main()
{

 printf(" %d \n", dist(-1.999999,0.0,0.0,0.0));
 return 0;
}
