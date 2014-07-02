// gcc m.c -Wall -lm
/*

for m=10.0000000000000000 : true = 0.0600000000000000; estim = 0.0551856465913334 ; ratio t/e =       1.0872392317
for m= 1.0000000000000000 : true = 0.0060000000000000; estim = 0.0005723179503263 ; ratio t/e =      10.4836830587 
for m= 0.1000000000000000 : true = 0.0006000000000000; estim = 0.0000057544987008 ; ratio t/e =     104.2662499727 
for m= 0.0100000000000000 : true = 0.0000600000000000; estim = 0.0000000575922671 ; ratio t/e =    1041.8065314555
for m= 0.0010000000000000 : true = 0.0000060000000000; estim = 0.0000000005767059 ; ratio t/e =   10403.9163074084  
for m= 0.0001000000000000 : true = 0.0000006000000000; estim = 0.0000000000057608 ; ratio t/e =  104152.4472835201  
for m= 0.0000100000000000 : true = 0.0000000600000000; estim = 0.0000000000000576 ; ratio t/e = 1041657.8145048500 

======================  +++++
for m=0.1000000000000000 : true = 0.0006000000000000; estim = 0.0000000000000000 ; ratio t/e =        inf 

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double pixelwidth = 0.006;
// DEM needs high bailout !!!!
double bailout = 1000.0;
unsigned int maxiter=1000;


int BadNumber(double d)
{ int r;

  switch (fpclassify(d))
 { 
  case FP_NAN:
  case FP_INFINITE:
  case FP_SUBNORMAL:
     r=1; // true
     break;
  case FP_ZERO:
  case FP_NORMAL:
   r= 0; // false 
  
 }
 
return r;

}



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
   double tmp;

   for (j = 1; j <= maxiter; j++)
    { tmp = x*xp - y*yp;
      if( BadNumber(tmp)) {printf ("for j = %d number x*xp - y*yp = %g ; number is bad  \n",j,tmp); break;}
       else nz = 2*tmp; 
      yp = 2*(x*yp + y*xp); 
      xp = nz; //zp = 2*z*zp;
      xp++; //zp = 2*z*zp + 1
      tmp =  x*x - y*y;
      if( BadNumber(tmp)) {printf ("for j = %d number x*x - y*y = %g ; number is bad  \n",j,tmp); break;}
         else nz = tmp + a; 
      y = 2*x*y + b; x = nz; //z = z*z + c;
      nz = x*x + y*y; 
      nzp = xp*xp + yp*yp;
      if (nzp > 1e40 ) return 0.0; // 
      if ( nz > bailout) break;
    }
   
   x = log(nz);
   if( BadNumber(x)) {printf ("x =log(nz) =  %g ; number is bad  \n",x); return x;} 
   x =x*x*nz/nzp;
   if( BadNumber(x)) {printf ("x =x*x*nz/nzp =  %g ; number is bad  \n",x); return x; }
   return fabs(log(x));
} //dist












int main()
{ 

  double td ; // true dist = pixelwidth*0.001;
  double ed; // estim dist 
  double m = 0.00000001;


 td = pixelwidth*m;
 ed=dist(-2.0-td,0.0,0.0,0.0);
 printf("for m=%.16f : true = %.16f; estim = %0.16f ; ratio t/e = %10.10f \n",m, td,ed , td/ed );
 
 
return 0;
}
