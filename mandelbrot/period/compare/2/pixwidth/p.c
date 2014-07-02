/*

gcc p.c -Wall -lm 
time ./a.out



numerical approximation of period of limit cycle 


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// THE REAL SLICE OF THE MANDELBROT SET 
long double CxMin = -1.4011552; //1890L; // The Feigenbaum Point = the limit of the period doubling cascade of bifurcations 
long double CxMax = 0.26L;  
long double Cx;
long double Cy = 0.0L;
long double PixelWidth ; // = (CxMax-CxMin)/10000.0L;
//long double precisionS ; //precisionS = PixelWidth / 100.0L;//= PixelWidth / 100.0L;
long double precisionJ = 1e-20; 
unsigned int periodJ, periodR;
long double Zp[2]; // periodic z points on dynamic plane


long double ER2 = 4.0L;
unsigned int jMax = 5000000; // iteration max = Max period 
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
   return (2*jMax+3); // (not escaping after 2*jMax = maybe periodic but period > jMax) or  
    // (maybe escaping but slow dynamics, so need more iterations then 2*jMax) 
}



// http://classes.yale.edu/Fractals/MandelSet/MandelScalings/CompDiam/CompDiam.html
unsigned int GivePeriodReal(long double Cx,long double Cy)
{
 long double Cx0= 0.25L; 
 long double Cx1= -0.75L; 
 long double Cx2= -1.25L; 
 long double Cx3= -1.368089448988708L; // numerical approximation = maybe wrong
 long double Cx4= -1.394040000725660L; // numerical approximation = maybe wrong
  
  if ( Cx1<Cx && Cx<Cx0 ) return 1;
  if ( Cx2<Cx && Cx<Cx1 ) return 2;
  if ( Cx3<Cx && Cx<Cx2 ) return 4; // numerical approximation = maybe wrong
  if ( Cx4<Cx && Cx<Cx3 ) return 8; // numerical approximation = maybe wrong
 return 0; // -1.36809742955000002314

}


// try to have the same number of the pixels = n
// inside each hyperbolic component of Mandelbrot set along real axis
// width of components 

long double GivePixelWidth(unsigned int period, unsigned int n)
{

  long double w ;
  unsigned int k;

 switch ( period )
 {  // A SCALING CONSTANT EQUAL TO UNITY IN 1D QUADRATIC MAPS M. ROMERA, G. PASTOR and F. MONTOYA
   case      0 : w=(CxMax-CxMin)/n;      break;
   case      1 : w=1.000000000000L/n;    break; // exact value
   case      2 : w=0.310700264133L/n;    break; // numerical approximation , maybe wrong 
   case      4 : w=0.070844843095L/n;    break; // w(2*p) = w(p)/4.0L  ; aproximated value of Feigenbaum constants
   case      8 : w=0.015397875272L/n;    break;
   case     16 : w=0.003307721510L/n;    break;
   case     32 : w=0.000708881730L/n;    break;
   case     64 : w=0.000151841994935L/n; break;
   case    128 : w=0.000032520887170L/n; break;
   case    256 : w=0.00000696502297L/n;  break;
   case    512 : w=0.000001491696694L/n; break;
   case   1024 : w=0.000000319475846L/n; break;
   case   2048 : w=0.000000068421948L/n; break;
   case   4096 : w=0.000000015L/n;       break;
   case   8192 : w=0.000000004L/n;       break;
   case  16384 : w=0.000000001L/n;       break;
   default : if (period == 2*jMax+3)  w=(CxMax-CxMin)/10.0L; // period not found or period > jMax
                else { k=period/16384; w = 0.000000001L; while (k>2) { w /=4.0L; k /=2;};  w /=n;} // dive 
 }

 return w;
}

int main()


{

PixelWidth = (CxMax-CxMin)/1000.0L;
precisionJ = PixelWidth/10000000.0L;
                    
// text file 
FILE * fp;  // result is saved to text file 
fp = fopen("data64_50.txt","w"); // create new file,give it a name and open it in binary mode  
fprintf(fp," periods of attracting orbits ( c points ) on real axis of parameter plane = real slice of the Mandelbrot set  \n");
fprintf(fp," from  Cmax = %.20Lf to Cmin = %.20Lf \n", CxMax, CxMin);
fprintf(fp," dC = CxMax-CxMin = %.20Lf \n", CxMax- CxMin);
fprintf(fp," non-inear scale with varied step = PixelWidth       \n");
fprintf(fp," precisionJ =  %.20Lf\n", sqrtl(precisionJ));
fprintf(fp,"  jMax = %u\n",  2*jMax);
fprintf(fp," \n\n\n");

// go along real axis from CxMin to CxMax using linear scale 
Cx = CxMax;
while (Cx>CxMin)
{ 
  // compute 
  periodR = GivePeriodReal(Cx,Cy);
  periodJ = GivePeriodJung(Cx, Cy, ER2, jMax, PixelWidth/10000000.0L, Zp);
  // check and save 
   if (periodJ == 2*jMax+3) 
      fprintf(fp," c = %.20Lf ; periodJ = %u ; PixelWidth = %.20LF Period not found : error !!! \n", Cx, periodJ, PixelWidth );
      else fprintf(fp," c = %.20Lf ; periodJ = %u ; PixelWidth = %.20LF \n", Cx, periodJ, PixelWidth );
  printf("c = %.20Lf ; period = %u \n",Cx, periodJ);  // info message
  // next c point 
  PixelWidth =GivePixelWidth( periodJ, 50);
  Cx -= PixelWidth;
}

 fclose(fp);
 printf(" result is saved to text file \n");




return 0;
}
 
