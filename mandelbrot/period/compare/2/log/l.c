/*

gcc l.c -Wall -lm 
time ./a.out



numerical approximation

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


long double ER2 = 4.0L;
unsigned int jMax = 500000; // iteration max = Max period 


// THE REAL SLICE OF THE MANDELBROT SET 
long double CxMin = -1.4011551890L; // The Feigenbaum Point = the limit of the period doubling cascade of bifurcations 
long double CxMax = -1.4L;
long double dC ; //= CxMax - CxMin;  
long double Cy = 0.0L;
// non linear scale of x axis 
long double lMax =  0.0L;  // logl(lMax) = 1.0
long double lMin =  -10.0L; // logl(lMin) = 0.0
long double dl ; // = lMax - lMin;
unsigned int iPixelsNumber = 10000; // number of points to check 
long double  step_l; //  = dl/iPixelsNumber;



                                   

long double GiveC( long double l)
{
  return CxMin + expl(l)*dC;
}





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






int main()


{

long double Cx;
long double l;


long double precisionJ = 1e-16; 
unsigned int  periodJ, periodR;
long double Zp[2]; // periodic z points on dynamic plane
unsigned int iMax = 2*jMax; // 1000000; // iteration max = Max period 
                    



// setup 
dC = CxMax - CxMin;  
dl = lMax - lMin;
step_l = dl/iPixelsNumber;



// text file 
FILE * fp;  // result is saved to text file 
fp = fopen("data10000.txt","w"); // create new file,give it a name and open it in binary mode  
fprintf(fp," periods of attracting orbits ( c points ) on real axis of parameter plane = real slice of the Mandelbrot set  \n");
fprintf(fp," from Cmin = %.20Lf to Cmax = %.20Lf \n", CxMin, CxMax);
fprintf(fp," dC = CxMax-CxMin = %.20Lf \n", dC);
fprintf(fp," Non-linear scale  \n");
fprintf(fp,"  precisionJ =  %.20Lf\n", sqrtl(precisionJ));
fprintf(fp," iMaxS = %u ; iMaxJ = %u\n", iMax, 2*jMax);
fprintf(fp," \n\n\n");




// go along real axis from CxMin to CxMax using nonlinear scale 
l = lMin;
Cx = GiveC(l);
printf("c = %.20Lf \n",Cx);
while (l<lMax)
{ 

  // compute 
  periodR = GivePeriodReal(Cx,Cy);
  periodJ = GivePeriodJung(Cx, Cy, ER2, jMax, precisionJ, Zp);
  // check and save 
  if (periodR>0)
    {
      if (periodJ==periodR ) // all periods are the same and real period is known 
         fprintf(fp," c = %.20Lf ; period = %u ; \n", Cx, periodJ );
         else fprintf(fp," c = %.20Lf ; period = %u ; periodJ = %u ; difference !!! \n", Cx, periodR, periodJ );
    }
    else // PeriodR==00
     {
       if (periodJ==0  ) 
         fprintf(fp," c = %.20Lf ; period = %u ; \n", Cx, periodJ );// all periods are the same and real period is known 
         else fprintf(fp," c = %.20Lf ;  periodJ = %u ;  \n", Cx, periodJ );
              
      }  

  
  // info message
  Cx = GiveC(l);
  printf("c = %.20Lf \n",Cx);
  // next  point 
  l += step_l;
}

   fclose(fp);
 printf(" result is saved to text file \n");
return 0;
}
 
