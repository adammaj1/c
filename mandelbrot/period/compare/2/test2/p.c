/*

gcc p.c -Wall -lm
time ./a.out



numerical approximation of period of limit cycle 


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



long double ER2 = 4.0L;
unsigned int jMax = 500000; 

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
unsigned int GivePeriodJung(long double cx, long double cy, long double ER2, unsigned int jMax, long double Zp[2]) 
{  //determine the period, then set Zp to periodic point.
   // bailout = ER2 = (EscapeRadius)^2
   unsigned int j;
  // unsigned int jMax = 500000; 
   long double x=0.0L;
   long double y=0.0L; // z 
   long double x0, y0; // z0 inside periodic orbit
   long double t; // temp
   long double precision = 1e-16;
   
   // iterate until z fall into periodic cycle ( = limit cycle) 
   for (j = 1; j <= jMax; j++)
   { 
     if (x*x + y*y <= ER2) 
       {t = x*x - y*y + cx; 
        y = 2*x*y + cy; 
        x = t;}
       else return 0; //escaping = definitely not periodic 
   } 
   // after jMax iterations z SHOULD BE inside periodic orbit 
   x0 = x; y0 = y; // z = z0

   // find period 
   for (j = 1; j <= jMax; j++)
   {  
      if (x*x + y*y <= ER2) 
        {t = x*x - y*y + cx; 
        y = 2*x*y + cy; 
        x = t;}
        else return 0; // escaping = definitely not periodic
      if ( (x - x0)*(x - x0) + (y - y0)*(y - y0) < precision)
      //if ( fabsl(x - x0)<precision && fabsl(y - y0)< precision) // periodic 
      {   Zp[0] = x; 
          Zp[1] = y; 
         return j;  // period = j 
      }
   }
   return (2*jMax+3); // not escaping after 2*jMax = maybe periodic 
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
unsigned int GivePeriodS(long double Cx,long double Cy, unsigned int iMax, long double precision, long double Zp[2])
{  
 
 
  long double Zx2, Zy2, /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
         ZPrevieousX,ZPrevieousY,
         ZNextX,ZNextY;
 
     unsigned int i; 
     unsigned int  period = iMax+3; // not periodic or period > iMax

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
long double CxMax = 0.26L;  
long double Cx;
long double Cy = 0.0L;
long double PixelWidth = (CxMax-CxMin)/10000.0L;
long double precisionS = PixelWidth / 100.0L;
long double precisionJ = sqrt(1e-16); 
unsigned int periodS, periodJ;
long double Zp[2]; // periodic z points on dynamic plane
unsigned int iMax = 1000000; // iteration max = Max period 
                    
// text file 
FILE * fp;  // result is saved to text file 
fp = fopen("data2p100.txt","w"); // create new file,give it a name and open it in binary mode  
fprintf(fp," periods of attracting orbits ( c points ) on real axis of parameter plane = real slice of the Mandelbrot set  \n");
fprintf(fp," from Cmin = %.20Lf to Cmax = %.20Lf \n", CxMin, CxMax);
fprintf(fp," dC = CxMax-CxMin = %.20Lf \n", CxMax- CxMin);
fprintf(fp," PixelWidth       = %.20Lf \n", PixelWidth);
fprintf(fp," precisionS        = %.20Lf ; precisionJ =  %.20Lf\n", precisionS, precisionJ);
fprintf(fp," iMaxS = %u ; iMaxJ = %u\n", iMax, 2*jMax);
fprintf(fp," \n\n\n");

// go along real axis from CxMin to CxMax
Cx = CxMin;
while (Cx<CxMax)
{

  periodS = GivePeriodS(Cx, Cy, iMax, precisionS, Zp);
  periodJ = GivePeriodJung(Cx, Cy, ER2, jMax, Zp);
  if (periodJ==periodS) 
         fprintf(fp," c = %.20Lf ; periodS = %u ; periodJ = %u \n", Cx, periodS, periodJ );
         else fprintf(fp," c = %.20Lf ; periodS = %u ; periodJ = %u ; difference !!! \n", Cx, periodS, periodJ );
  printf("c = %.20Lf \n",Cx);
  Cx += PixelWidth;
}

 fclose(fp);
 printf(" result is saved to text file \n");
return 0;
}
 
