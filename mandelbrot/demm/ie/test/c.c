#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


double EscapeRadius = 33.0;  /* radius of circle around origin; its complement is a target set for escaping points */
double ER2 ; //= (EscapeRadius*EscapeRadius)

unsigned int GivePeriodJung( long double cx, long double cy,  long double ER2,  long double Zp[2]) 
{  //determine the period, then set Zp to periodic point.
   // bailout = ER2 = (EscapeRadius)^2
   unsigned int j;
   unsigned int jMax = 50000; 
   long double x=0.0L;
   long double y=0.0L; // z 
   long double x0, y0; // z0 inside periodic orbit
   long double t; // temp
   long double precision = 1e-16;
   
   // iterate until z fall into periodic cycle 
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
      if ( (x - x0)*(x - x0) + (y - y0)*(y - y0) < 1e-16)
      //if ( fabsl(x - x0)<precision && fabsl(y - y0)< precision) // periodic 
      {   Zp[0] = x; 
          Zp[1] = y; 
         return j;  // period = j 
      }
   }
   return 10000; // not escaping after 2*jMax = maybe periodic 
}



/*-----------------------------*/
int SameComplexValue(double Z1x,double Z1y,double Z2x,double Z2y, double precision)
{
    if (fabs(Z1x-Z2x)<precision && fabs(Z1y-Z2y)<precision) 
       return 1; /* true */
       else return 0; /* false */
    }
 
/*-------------------------------*/
// this function is based on program:
// Program MANCHAOS.BAS  
// http://sprott.physics.wisc.edu/chaos/manchaos.bas
// (c) 1997 by J. C. Sprott 
//
unsigned int GivePeriodS(long double Cx, long double Cy, int Iteration_Max, long double precision, long double Zp[2])
{  
 
 
  long double Zx2, Zy2, /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
         ZPrevieousX,ZPrevieousY,
         ZNextX,ZNextY;
 
     int Iteration,
      I; 
     int  period = 0; // not periodic or period > Iteration_Max
  long double orbit[Iteration_Max+1][2]; /* array elements are numbered from 0 to length-1 */    
 
  Zp[0] = 0.0;
  Zp[1] = 0.0; 

  /* starting point is critical point  */
   ZPrevieousX=0.0;
   ZPrevieousY=0.0;
   orbit[0][0]=0.0;
   orbit[0][1]=0.0;  
   Zx2=ZPrevieousX*ZPrevieousX;
   Zy2=ZPrevieousY*ZPrevieousY;
   /* iterate and save points for analysis */
   for (Iteration=1;Iteration<Iteration_Max+1 ;Iteration++)
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
            orbit[Iteration][0]=ZNextX;
            orbit[Iteration][1]=ZNextY;   
 
        }; 
    /* here iteration=IterationMax+1 but last element of orbit has number IterationMax */    
     for(I=Iteration_Max-1;I>0;I--) 
      if (SameComplexValue(orbit[Iteration_Max][0],orbit[Iteration_Max][1],orbit[I][0],orbit[I][1],precision))
        { 
          Zp[0] = orbit[Iteration_Max][0];
          Zp[1] = orbit[Iteration_Max][1]; 
          period = Iteration_Max-I;
          break;
        }
   
  
  return period ; 
}
 
 /**
   * Calculates a lower bound for the distance between C 
   * and the closest point in the border of M.
   * @param cycle double[][] - C's periodic cycle. = double Zx0, double Zy0
   * @return double - The interior distance estimation.
    http://www.moleculardensity.net/buddhabrot/appendix/2
    Interior and exterior distance bounds for the Mandelbrot set by Albert Lobo Cusidó 
   */
   double getInteriorLowerBound(int period, double Zp[2]) {
    // Real and imaginary components for complex numbers D1, D2, D3, and D4;
    // and temporary variables to store the complex numbers in recursion.

    double Zr, Zi, D1r, D1i, D2r, D2i, D3r, D3i, D4r, D4i;
    double D1rT, D1iT, D2rT, D2iT, D3rT, D3iT, D4rT, D4iT;
    int i;
    double distance;

    // Initial values:  D1 = 1;  D2 = 0;  D3 = 0;  D4 = 0;
    D1r = 1;
    D1i = D2r = D2i = D3r = D3i = D4r = D4i = 0;

    // Start iterating.
    for (i = 0; i < period; i++) {
      // No need to iterate Z; values are already in 'cycle'.

      Zr = Zp[0];
      Zi = Zp[1];

      // D1 = 2 * Z * D1;
      D1rT = 2 * (Zr * D1r - Zi * D1i);
      D1iT = 2 * (Zi * D1r + Zr * D1i);

      // D2 = 2 * Z * D2 + 1;
      D2rT = 2 * (Zr * D2r - Zi * D2i) + 1;
      D2iT = 2 * (Zi * D2r + Zr * D2i);

      // D3 = 2 * (D1^2 + Z * D3);
      D3rT = 2 * ((Zr * D3r - Zi * D3i) + (D1r * D1r - D1i * D1i));
      D3iT = 2 * ((Zr * D3i + Zi * D3r) + (2 * D1r * D1i));

      // D4 = 2 * (D1 * D2 + Z * D4);
      D4rT = 2 * ((Zr * D4r - Zi * D4i) + (D1r * D2r - D1i * D2i));
      D4iT = 2 * ((Zr * D4i + Zi * D4r) + (D1r * D2i + D1i * D2r));

      // Update variables.

      D1r = D1rT;  D1i = D1iT;
      D2r = D2rT;  D2i = D2iT;
      D3r = D3rT;  D3i = D3iT;
      D4r = D4rT;  D4i = D4iT;
    }

    // A = 1 - |D1|^2;
    double A = 1 - (D1r * D1r + D1i * D1i);
    // B = |D4 + D3 * (D2 / (1 - D1))|;
    double B = (1 - D1r) * (1 - D1r) + D1i * D1i;
    D1rT = (D2r * (1 - D1r) - D2i * D1i) / B;
    D1iT = (D2i * (1 - D1r) + D2r * D1i) / B;
    D2rT = D4r + (D3r * D1rT - D3i * D1iT);
    D2iT = D4i + (D3i * D1rT + D3r * D1iT);
    B    = sqrt(D2rT * D2rT + D2iT * D2iT);

    // Return lower bound -that is, 1/4 the estimated bound.
    distance =  log(A / (4.0 * B));
    

    
    return distance ;
  }

/*
 estimates distance from point c to nearest point in Julia  set 
 for Fc(z)= z*z + c
 z(n+1) = Fc(zn)  
 this function is based on function  mndlbrot::dist  from  mndlbrot.cpp
 from program mandel by Wolf Jung (GNU GPL )
 http://www.mndynamics.com/indexp.html 

Hyunsuk Kim  : 
For Julia sets, z is the variable and c is a constant. Therefore df[n+1](z)/dz = 2*f[n]*f'[n] -- you don't add 1.

For the Mandelbrot set on the parameter plane, you start at z=0 and c becomes the variable. df[n+1](c)/dc = 2*f[n]*f'[n] + 1. 
http://iquilezles.org/www/articles/distancefractals/distancefractals.htm

double EscapeRadius = 33.0; 

 */

// boolean Escape time and DEM/M in one loop  
double GiveDistance( double Cx, double Cy, int iMax, double DistanceMax)
{ 

  // C = Cx + Cy* I = point of parameter c-plane

  int i; /* iteration number */

  double Zx, Zy; /* Z = Zx + Zy*I  point of dynamical plane */
  double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
  double temp;
  double absZ2;

  // http://en.wikipedia.org/wiki/Complex_quadratic_polynomial
  // first derivative  of fc(zcr) with respect to c =  dZ = dfc(zcr)/dc = 2*Z*dZ  = dZx + dZy*I 
  double dZx = 0.0;
  double dZy = 0.0; 
  double absdZ2; // = abs(dZ)* abs(dZ) = dZx*dZx + dZy*dZy

  double distance;
  

  /* initial value of orbit  = critical point zcr = 0 */
  Zx=0.0; 
  Zy=0.0;
  //
  Zx2 = Zx*Zx;
  Zy2 = Zy*Zy;
  absZ2= Zx*Zx + Zy*Zy;
  absdZ2= dZx*dZx + dZy*dZy;
  // iteration of critical point z= 0 on the dynamical z-plane       
  for (i=0; i<iMax; i++)
    { // check if not escaping : abs(z)>ER
      if (absZ2  > ER2 )  break ;  // exterior when escapes
      //if (absdZ2 > 1e60) { i=iMax;  } //  interior when derivative explodes 
     

       // in the loop, the derivative should be calculated before the new z
      /* first derivative   zp = 2*z*zp  = xp + yp*i; */
      temp = 2*(Zx*dZx - Zy*dZy) + 1.0 ; /*  */ 
      dZy = 2*(Zx*dZy + Zy*dZx); 
      dZx = temp;

      // z = fc(z) = z*z + c 
      Zy=2*Zx*Zy + Cy;
      Zx=Zx2-Zy2 +Cx;
      
     

      // abs 
      Zx2 = Zx*Zx;      
      Zy2 = Zy*Zy;
      absZ2= Zx2 + Zy2; // nz = x*x + y*y;
      absdZ2= dZx*dZx + dZy*dZy; // 
    };

 // compute distance 
  if (i<iMax) // exterior
   { 
     
     distance = sqrt(absZ2/absdZ2)*log(absZ2); //
     distance = pow( 4.0*distance, 0.25 );
                 
    }
   else distance=0.0; // interior 
   // if (nz < bailout) return 1; //still not escaping after iteration , rare
  // if (absdZ2 < absZ2) color= iExterior;  //includes escaping through 0 // som eiterior points are coded as exterior = error
    
  return distance; 
} //0.000007
 

int main()


{
// A real period 5 hyperbolic component whose center is 
// approximately located at c=−1.985424253. 
// external parameter angles :  15/31 , 16/33
// EXTERNAL RAYS AND THE REAL SLICE OF THE MANDELBROT SET by SAEED ZAKERI
// 
long double CxMin = -1.401155L; 
long double CxMax = 0.26L;  
long double Cx;
long double Cy = 0.0L;
long double PixelWidth = (CxMax-CxMin)/10000.0L;
long double precision = PixelWidth / 1.0000000000L; 
unsigned int periodS, periodJ;
long double Zp[2]; // periodic z points on dynamic plane
unsigned int iMax = 100000; // iteration max = Max period 
                    
// text file 
FILE * fp;  // result is saved to text file 
fp = fopen("data2.txt","w"); /*create new file,give it a name and open it in binary mode  */
fprintf(fp," periods of attracting orbits ( c points ) on real axis of parameter plane = real slice of the Mandelbrot set  \n");
fprintf(fp," from Cmin = %.20Lf to Cmax = %.20Lf \n", CxMin, CxMax);
fprintf(fp," dC = CxMax-CxMin = %.20Lf \n", CxMax- CxMin);
fprintf(fp," PixelWidth       = %.20Lf \n", PixelWidth);
fprintf(fp," precision        = %.20Lf \n", precision);
fprintf(fp," iMax = %u \n", iMax);
fprintf(fp," \n\n\n");

// go along real axis from CxMin to CxMax
Cx = CxMin;
while (Cx<CxMax)
{

  periodS = GivePeriodS(Cx, Cy, iMax, precision, Zp);
  periodJ = GivePeriodJung(Cx, Cy, ER2, Zp);
  fprintf(fp," c = %.20Lf ; periodS = %u ; periodJ = %u \n", Cx, periodS, periodJ );
  Cx += PixelWidth;
}

 fclose(fp);
 printf(" result is saved to text file \n");
return 0;
}
 
