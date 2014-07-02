/*

gcc p.c -Wall -lm 
time ./a.out



numerical approximation  of limit cycle's period  
along real slice of Mandelbrot set  

Adam Majewski


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// part of THE REAL SLICE OF THE MANDELBROT SET where period doubling cascade is  
long double CxMin ; //= -1.4011552; // 1890L; // > The Feigenbaum Point = the limit of the period doubling cascade of bifurcations 
long double CxMax ; //= 0.26L;  
long double Cx;
long double Cy = 0.0L; // constant value 
long double PixelWidth ; // = (CxMax-CxMin)/10000.0L;

 
long double precision2 = 1e-20; 
unsigned int periodJ, periodR;
long double Zp[2]; // periodic z points on dynamic plane


long double ER2 = 4.0L;
unsigned int ixMax = 200; // number of points to check 
long double dxMin = 0.0L; // 
long double dxMax = 30.0L; // 
long double dxStep ; //= (dxMax-dxMin)/((long double)ixMax);
unsigned int jMax = 5000000; // iteration max = Max period 
unsigned int iNoPeriod;
//unsigned int iMax ; //= 2*jMax; // 1000000; // iteration max = Max period 

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

code with small changes 

*/
unsigned int GivePeriodJung(long double cx, long double cy, long double ER2, unsigned int jMax, long double precision2, long double Zp[2]) 
{  //determine the period, then set Zp to periodic point.
   // bailout = ER2 = (EscapeRadius)^2
   unsigned int j;
  // unsigned int jMax = 500000; 
   long double x=0.0L;
   long double y=0.0L; // z 
   long double x0, y0; // z0 inside periodic orbit
   long double t; // temp
   //long double precision = 1e-16;
   
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

   // find a period 
   for (j = 1; j <= jMax; j++)
   {  
      if (x*x + y*y <= ER2) 
        {t = x*x - y*y + cx; 
        y = 2*x*y + cy; 
        x = t;}
        else return 0; // escaping = definitely not periodic
      
     if ( (x - x0)*(x - x0) + (y - y0)*(y - y0) < precision2) // periodic 
      {   Zp[0] = x; 
          Zp[1] = y; 
         return j;  // period = j 
      }
   }
   return (iNoPeriod); // (not escaping after 2*jMax = maybe periodic but period > jMax) or  
    // (maybe escaping but slow dynamics, so need more iterations then 2*jMax) 
}



long double Givedx(unsigned int i)
{
 return dxMin + i*dxStep;
}


// approiximated function using roots on the real axis
long double GiveCx(long double x)
{ // http://zunzun.com
  
  long double a = 0.53226927610784935L;
  long double b = 0.65410208763684241L;
  long double c = -1.4312869957125389L;
  long double d = 0.84710834303177074L;
  return (c*atanl(expl((x-a)/b)) + d);
}




int main()
{


iNoPeriod = 2*jMax+3;
dxStep = (dxMax-dxMin)/((long double)ixMax);
CxMin = GiveCx(dxMax);
CxMax = GiveCx(dxMin);

int ix;
long double dx;
long double oldCx =CxMax + 0.003L;
                    
// text file 
FILE * fp;  // result is saved to text file 
fp = fopen("data30.txt","w"); // create new file,give it a name and open it in binary mode  
fprintf(fp," periods of attracting orbits ( c points ) on real axis of parameter plane = real slice of the Mandelbrot set  \n");
fprintf(fp," from  Cmax = %.20Lf to Cmin = %.20Lf \n", CxMax, CxMin);
fprintf(fp," dC = CxMax-CxMin = %.20Lf \n", CxMax- CxMin);
fprintf(fp," non-inear scale with varied step = PixelWidth       \n");
fprintf(fp," precision =  %.20Lf\n", sqrtl(precision2));
fprintf(fp,"  jMax = %u\n",  2*jMax);
fprintf(fp," \n\n\n");

// go along real axis from CxMin to CxMax using nonlinear scale 
for (ix=0; ix<ixMax; ix++)

{ 
  // compute
  dx = Givedx(ix);
  Cx= GiveCx(dx);
  PixelWidth = oldCx - Cx;
  oldCx = Cx; 
  //periodR = GivePeriodReal(Cx,Cy);
  periodJ = GivePeriodJung(Cx, Cy, ER2, jMax, precision2, Zp);
  // check and save 
   if (periodJ == iNoPeriod) 
      fprintf(fp," c = %.20Lf ; periodJ = %u ; PixelWidth = %.20LF Period not found : not escaping or not periodic after %d !!! \n", Cx, periodJ, PixelWidth , 2*jMax);
      else fprintf(fp," c = %.20Lf ; periodJ = %u ; PixelWidth = %.20LF \n", Cx, periodJ, PixelWidth );
  printf("c = %.20Lf ; period = %u \n", Cx, periodJ);  // info message
  
 
}

 fclose(fp);
 printf(" result is saved to text file \n");





return 0;
}
 
