/*

gcc pc -Wall -lm
time ./a.out






*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



long double ER2 = 4.0L;

// mndynamics::period(double &a, double &b, int cycle)
// mndynamo.cpp  by Wolf Jung (C) 2007-2014
// part of Mandel 5.10 which is free software; you can
//   redistribute and / or modify them under the terms of the GNU General
//   Public License as published by the Free Software Foundation; either
//   version 3, or (at your option) any later version. In short: there is
//   no warranty of any kind; you must redistribute the source code as well.
/*

void mndlbrot::f(double a, double b, double &x, double &y) const
{ double u = x*x - y*y + a; y = 2*x*y + b; x = u; }
*/
int GivePeriodJung(double &a, double &b, int cycle) // = 0
{  //determine the period, if cycle then set a, b to periodic point.
   uint j; 
   double x, y, x0, y0; 
   
   critical(a, b, x, y);
   //
   for (j = 1; j <= 1000; j++)
   { if (x*x + y*y <= bailout) f(a, b, x, y); else return 0; }
   x0 = x; y0 = y; // after jMax iterations z should be inside periodic orbit 
   //
   for (j = 1; j <= 1024; j++)
   {  if (x*x + y*y <= bailout) f(a, b, x, y); else return 0;
      if ( (x - x0)*(x - x0) + (y - y0)*(y - y0) < 1e-16)
      {  if (cycle) { a = x; b = y; }
         return j;
      }
   }
   return 10000; // 
}

int SameComplexValue(long double Z1x,long double Z1y,long double Z2x,long double Z2y, long double precision)
{
    if (fabsl(Z1x-Z2x)<precision && fabs(Z1y-Z2y)<precision) 
       return 1; /* true */
       else return 0; /* false */
    }
 
/*-------------------------------*/
// this function is based on program:
// Program MANCHAOS.BAS  
// http://sprott.physics.wisc.edu/chaos/manchaos.bas
// (c) 1997 by J. C. Sprott 
//
int GivePeriodS(long double Cx,long double Cy, int iMax, long double precision, long double Zp[2])
{  
 
 
  long double Zx2, Zy2, /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
         ZPrevieousX,ZPrevieousY,
         ZNextX,ZNextY;
 
     int i; 
     int  period = 0; // not periodic or period > iMax

     /* dynamic 1D arrays for  x, y of z points   */
    long double *OrbitX; // zx
    long double *OrbitY;  // zy 
     int iLength = iMax; // length of arrays ;  array elements are numbered from 0 to iMax-1 
  //  creates dynamic arrays and checks if it was done properly
  OrbitX = malloc( iLength * sizeof(long double) );
  OrbitY = malloc( iLength * sizeof(long double) );
  if (OrbitX == NULL || OrbitY ==NULL)
    {
      printf("Could not allocate memory \n");
      return 1; // error
    }

 
  Zp[0] = 0.0;
  Zp[1] = 0.0; 

  /* starting point is critical point  */
   ZPrevieousX=0.0;
   ZPrevieousY=0.0;
   OrbitX[0] =0.0;
   OrbitY[0] =0.0;  
   Zx2=ZPrevieousX*ZPrevieousX;
   Zy2=ZPrevieousY*ZPrevieousY;

   /* iterate and save points to the array */
   for (i=0;i<iMax ;i++)
        {
            ZNextY=2*ZPrevieousX*ZPrevieousY + Cy;
            ZNextX=Zx2-Zy2 +Cx;
            Zx2=ZNextX*ZNextX;
            Zy2=ZNextY*ZNextY;
            if ((Zx2+Zy2)>ER2) return 0; /* basin of atraction to infinity */
            //if (SameComplexValue(ZPrevieousX,ZPrevieousY,ZNextX,ZNextY,precision))
            //   return 1; /* fixed point , period =1 */
            ZPrevieousX=ZNextX;
            ZPrevieousY=ZNextY;
            /* */
            OrbitX[i] = ZNextX;
            OrbitY[i] = ZNextY;   
 
        };
 
    /* find   */    
     for(i=iMax-2;i>0;i--) 
      if (SameComplexValue(OrbitX[iMax-1],OrbitY[iMax-1],OrbitX[i],OrbitY[i],precision))
        { 
          Zp[0] = OrbitX[i];
          Zp[1] = OrbitY[i]; 
          period = iMax-i-1; // compute period 
          break; // the loop 
        }
   
  // free memmory
  free(OrbitX);
  free(OrbitY);

  return period ; 
}

int main()


{
// A real period 5 hyperbolic component whose center is 
// approximately located at c=−1.985424253. 
// external parameter angles :  15/31 , 16/33
// EXTERNAL RAYS AND THE REAL SLICE OF THE MANDELBROT SET by SAEED ZAKERI
// 
long double CxMin = -1.9855056L; 
long double CxMax = -1.98538L;  
long double Cx;
long double Cy = 0.0L;
long double PixelWidth = (CxMax-CxMin)/10000.0L;
long double precision = PixelWidth / 10.0L; 
int periodS;
long double Zp[2]; // periodic z points on dynmaic plane
int iMax = 1000000;

// text file 
FILE * fp;  // result is saved to text file 
fp = fopen("data5.txt","w"); /*create new file,give it a name and open it in binary mode  */
fprintf(fp," periods of c points on real axis of parameter plane = real slice of the Mandelbrot set  \n");
fprintf(fp," from Cmin = %.20Lf to Cmax = %.20Lf \n", CxMin, CxMax);
fprintf(fp," dC = CxMax-CxMin = %.20Lf \n", CxMax- CxMin);
fprintf(fp," PixelWidth       = %.20Lf \n", PixelWidth);
fprintf(fp," precision        = %.20Lf \n", precision);
fprintf(fp," iMax = %d \n", iMax);
fprintf(fp," \n\n\n");

// go along real axis from CxMin to CxMax
Cx = CxMin;
while (Cx<CxMax)
{

  periodS = GivePeriodS(Cx, Cy, iMax, precision, Zp);
  fprintf(fp," c = %.20Lf ; periodS = %d \n", Cx, periodS );
  Cx += PixelWidth;
}

 fclose(fp);
 printf(" result is saved to text file \n");
return 0;
}
 
