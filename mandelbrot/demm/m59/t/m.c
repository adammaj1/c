// gcc m.c -Wall -lm
/*

true = 0.600000; estim = 0.8554400314182631 ; ratio t/e = 0.701393
true = 0.060000; estim = 0.0007582245426766 ; ratio t/e = 79.132231 
true = 0.006000; estim = 0.0000007318820593 ; ratio t/e = 8198.042189
true = 0.000600; estim = 0.0000000007291222 ; ratio t/e = 822907.323994  
true = 0.000060; estim = 0.0000000000007288 ; ratio t/e = 82322007.813895 
true = 0.000006; estim = 0.0000002008056765 ; ratio t/e = 29.879633  
true = 0.000001; estim = 0.0000035929010562 ; ratio t/e = 0.166996 

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double pixelwidth = 0.006;
// DEM needs high bailout !!!!
double bailout = 1000.0;
unsigned int maxiter=1000;

/*

if (nzp < 0.1 || nz < 1) return (255 << 16) | (255 << 8) | 0; // boundary
  x = log(nz);
  if (x*x*nz/nzp < step_x/1000) return (255 << 16) | (0 << 8) | 0; // boundary - most points
     else return (255 << 16) | (255 << 8) | 255; // exterior
*/

// mndlbrot.cpp from Mandel 5.9 by Wolf Jung (C) 2007-2013
// c=a+b*i z=x+y*i
double dist(double a, double b, double x, double y)
{  

   unsigned int j; 
   double xp = 1; //zp = 1
   double yp = 0;
   double nz=0.0;
   double  nzp=0.0; 

   for (j = 1; j <= maxiter; j++)
    { 

      nz = 2*(x*xp - y*yp); 
      yp = 2*(x*yp + y*xp); 
      xp = nz; //zp = 2*z*zp;
      xp++; //zp = 2*z*zp + 1
      nz = x*x - y*y + a; 
      y = 2*x*y + b; x = nz; //z = z*z + c;
      nz = x*x + y*y; 
      nzp = xp*xp + yp*yp;
      if (nzp > 1e40 ) return 0.0;
      if ( nz > bailout) break;
    }
   //if (nz < bailout) return 1; //not escaping, rare
   //if (nzp < nz) return 10; //includes escaping through 0
   // d = sqrt( dot(z,z)/dot(dz,dz) )*log(dot(z,z));
   x = log(nz); 
   x =x*x*nz/nzp;
   return x;
} //dist


int main()
{ double d = pixelwidth*0.001;

 printf(" true = %f; estim = %.16f ; ratio t/e = %f  \n",d, dist(0.25+d,0.0,0.0,0.0), d/dist(0.25+d,0.0,0.0,0.0));
 return 0;
}
