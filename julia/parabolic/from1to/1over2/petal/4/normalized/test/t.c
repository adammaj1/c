/*


gcc t.c -lm -Wall
./a.out

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>


const int iMax =20000;
//int iMaxN;

/* fc(z) = z*z + c */
 int iPeriodChild ; //// Period of secondary component joined by root point with the parent component 
int iPeriodParent = 1; // main cardioid of Mandelbrot set
unsigned int denominator;/* denominator of internal angle */
double InternalAngle;
// c parameter of the function fc(z)
double Cx;
double Cy; 
double complex c;
// parabolic fixed point
double complex ZA;  /* atractor ZA = ZAx + ZAy*i */
double ZAx,ZAy;
// critical point
double complex ZCr = 0.0; 
double ZCrx = 0.0;
double ZCry = 0.0; 



double AR ; //=0.05  /* PixelWidth*1.5   radius of circle around attractor ZA = target set for attracting points */
double  AR2 ; //AR*AR
double ARn; // normalized see : GiveNormalizedParameters ; normalization : phi(critical point ) = 0
double ARn2;

double EscapeRadius=2.0; /* radius of circle around origin; its complement is a target set for escaping points */
double ER2; //=EscapeRadius*EscapeRadius;
//#define alfa (1-sqrt(1-4*Cx))/2 /* attracting or parabolic fixed point z = alfa */
//#define beta (1+sqrt(1-4*Cx))/2 /* repelling or parabolic fixed point z = beta */




// integer ( screen) coordinate of virtual 2D array 
// Indexes of array starts from 0 not 1 
//  unsigned int ixMin = 0; // Indexes of array starts from 0 not 1
unsigned int iXmax ; //
static unsigned int iWidth = 2000; // horizontal dimension of array
//static unsigned int iyMin = 0; // Indexes of array starts from 0 not 1
unsigned int iYmax ; //
static unsigned int iHeight = 1000; //  odd number !!!!!! = (iyMax -iyMin + 1) = iyAboveAxisLength + iyBelowAxisLength +1
// The size of array has to be a positive constant integer 
unsigned int  iLength;  


/* world ( double) coordinate = parameter plane*/
const double ZxMin=-2.0;
const double ZxMax=2.0;
const double ZyMin=-1.0;
const double ZyMax=1.0;
double PixelWidth; //=(ZxMax-ZxMin)/iXmax;
double PixelHeight; //=(ZyMax-ZyMin)/iYmax;


/* find c in component of Mandelbrot set 
 
   uses code by Wolf Jung from program Mandel
   see function mndlbrot::bifurcate from mandelbrot.cpp
   http://www.mndynamics.com/indexp.html
 
*/
double complex GiveC(double InternalAngleInTurns, double InternalRadius, unsigned int iPeriodOfTheParent)
{
  //0 <= InternalRay<= 1
  //0 <= InternalAngleInTurns <=1
  double t = InternalAngleInTurns *2*M_PI; // from turns to radians
  double R2 = InternalRadius * InternalRadius;
  //double Cx, Cy; /* C = Cx+Cy*i */
  switch ( iPeriodOfTheParent ) // of component 
    {
    case 1: // main cardioid
      Cx = (cos(t)*InternalRadius)/2-(cos(2*t)*R2)/4; 
      Cy = (sin(t)*InternalRadius)/2-(sin(2*t)*R2)/4; 
      break;
    case 2: // only one component 
      Cx = InternalRadius * 0.25*cos(t) - 1.0;
      Cy = InternalRadius * 0.25*sin(t); 
      break;
      // for each iPeriodOfTheParent  there are 2^(iPeriodOfTheParent-1) roots. 
    default: // higher periods : not works, to do
      Cx = 0.0;
      Cy = 0.0; 
      break; }
 
  return Cx + Cy*I;
}
 

/*
 
  http://en.wikipedia.org/wiki/Periodic_points_of_complex_quadratic_mappings
  z^2 + c = z
  z^2 - z + c = 0
  ax^2 +bx + c =0 // ge3neral for  of quadratic equation
  so :
  a=1
  b =-1
  c = c
  so :
 
  The discriminant is the  d=b^2- 4ac 
 
  d=1-4c = dx+dy*i
  r(d)=sqrt(dx^2 + dy^2)
  sqrt(d) = sqrt((r+dx)/2)+-sqrt((r-dx)/2)*i = sx +- sy*i
 
  x1=(1+sqrt(d))/2 = beta = (1+sx+sy*i)/2
 
  x2=(1-sqrt(d))/2 = alfa = (1-sx -sy*i)/2
 
  alfa : attracting when c is in main cardioid of Mandelbrot set, then it is in interior of Filled-in Julia set, 
  it means belongs to Fatou set ( strictly to basin of attraction of finite fixed point )
 
*/
// uses global variables : 
//  ax, ay (output = alfa(c)) 
double complex GiveAlfaFixedPoint(double complex c)
{
  double dx, dy; //The discriminant is the  d=b^2- 4ac = dx+dy*i
  double r; // r(d)=sqrt(dx^2 + dy^2)
  double sx, sy; // s = sqrt(d) = sqrt((r+dx)/2)+-sqrt((r-dx)/2)*i = sx + sy*i
  double ax, ay;
 
  // d=1-4c = dx+dy*i
  dx = 1 - 4*creal(c);
  dy = -4 * cimag(c);
  // r(d)=sqrt(dx^2 + dy^2)
  r = sqrt(dx*dx + dy*dy);
  //sqrt(d) = s =sx +sy*i
  sx = sqrt((r+dx)/2);
  sy = sqrt((r-dx)/2);
  // alfa = ax +ay*i = (1-sqrt(d))/2 = (1-sx + sy*i)/2
  ax = 0.5 - sx/2.0;
  ay =  sy/2.0;
 
 
  return ax+ay*I;
}


/*
 normalization : 
   phi(critical point ) = 0
    Attraction time(zcr) = k
   attraction time of (zcr+pixelwidth) = k+1
  so :
fix iMax
 compute k
 compute distance between (zcr+pixelwidth) and alpha after k iterations 
 change AR : d( zcr, alfa, k) < newAR < d( zcr+pixelWidth,alfa, k)

then zcr will be on boundary of level sets 

 here z_cr = 0.0
 changes iterationMax and AR 
*/
double GiveNormalizedParameters()
{
 
 int i,p;
 int iMaxN;
 
// initiall point of iteration is a critical point of f
 double Zx=0.0;
 double  Zy=0.0;
 double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
 /* distance from z to Alpha  after iMax iterations */
 double dcr2, dcr;
 double  d1;
 
 // compute dcr
 Zx = ZCrx; // critical point z=0
 Zy = ZCry;
 Zx2=Zx*Zx;
 Zy2=Zy*Zy;
  
 for ( i=0; i<iMax; ++i)
  {   //  forward iteration of z under fc(z) = z^2 + c
      for(p=0; p< iPeriodChild ;++p) // iMax = period !!!!
       { 
        // not needed here if (Zx2+Zy2 > ER2) return iExterior; // bailout test 
        Zy=2*Zx*Zy + Cy;
        Zx=Zx2-Zy2 +Cx;
        Zx2=Zx*Zx;
        Zy2=Zy*Zy;
        
       }
  dcr2 = (Zx-ZAx)*(Zx- ZAx) + (Zy-ZAy)*(Zy-ZAy); // 
  if (dcr2<AR2)  {iMaxN = i; break;}

  }
  
 
  

 // compute d1 // next pixel from critical point 
 Zx = ZCrx; 
 Zy = ZCry+PixelHeight; // move up or down , because next right/left  point has the same time  
 Zx2=Zx*Zx;
 Zy2=Zy*Zy;
 for ( i=0; i<iMaxN+1; ++i)
  {   //  forward iteration of z under fc(z) = z^2 + c
      for(p=0; p< iPeriodChild ;++p) // iMax = period !!!!
       { 
        // not needed here if (Zx2+Zy2 > ER2) return iExterior; // bailout test 
        Zy=2*Zx*Zy + Cy;
        Zx=Zx2-Zy2 +Cx;
        Zx2=Zx*Zx;
        Zy2=Zy*Zy;
         
      }
            
      
  }




  d1 = sqrt((Zx-ZAx)*(Zx- ZAx) + (Zy-ZAy)*(Zy-ZAy)); 
  dcr = sqrt(dcr2);

  if (dcr<AR && AR<d1 && d1<dcr ) 
   { printf("  good parameters, nothing to do : d1 = %.20f < dcr = %.20f < AR = %.20f \n", d1,dcr, AR);
     ARn = AR;
  }
   else 
{   

     // now : dcr<d1<AR
     // find ARn = normalised AR  
    ARn = dcr - (dcr-d1)*0.5;
    // now dcr<AR<d1
    printf(" normalised parameters : d1 = %.20f < dcr = %.20f < ARn = %.20f ;  AR = %.20f  \n",d1,  dcr, ARn, AR);
 } 
 return ARn; // 
 
}





int GiveAttractiveTime(double complex z)
{
 
 int i,p;
 //int iMaxN;
 
// initiall point of iteration 
 double Zx=creal(z);
 double Zy=cimag(z);
 double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
 /* distance from z to Alpha  after iMax iterations */
 double d2;
  
 // compute dcr
 Zx2=Zx*Zx;
 Zy2=Zy*Zy;
  
 for ( i=0; i<iMax; ++i)
  {   //  forward iteration of z under fc(z) = z^2 + c
      for(p=0; p< iPeriodChild ;++p) // iMax = period !!!!
       { 
        // not needed here if (Zx2+Zy2 > ER2) return iExterior; // bailout test 
        Zy=2*Zx*Zy + Cy;
        Zx=Zx2-Zy2 +Cx;
        Zx2=Zx*Zx;
        Zy2=Zy*Zy;
        
       }
  d2 = (Zx-ZAx)*(Zx- ZAx) + (Zy-ZAy)*(Zy-ZAy); // 
  if (d2<ARn2)  break;

  }
 printf(" z = %f ; %f ; d = %0.16f ", creal(z), cimag(z), sqrt(d2));
  
 return i; // 
 
}








int setup()
{

iLength = iWidth*iHeight;/* length of array in bytes = number of bytes = number of pixels of image * number of bytes of color */
 


  PixelWidth=(ZxMax-ZxMin)/ (iWidth-1);
  PixelHeight=(ZyMax-ZyMin)/(iHeight-1);
  
  
 

  denominator = iPeriodChild;
  InternalAngle = 1.0/((double) denominator);
  // find the c parameter with specyfic features
  c = GiveC(InternalAngle, 1.0, iPeriodParent) ; // internal radius= 1.0 gives root point = parabolic parameter   
  Cx=creal(c);
  Cy=cimag(c);
 
  
 ZA = GiveAlfaFixedPoint( c);
 ZAx = creal(ZA);
 ZAy = cimag(ZA);

 ER2=EscapeRadius*EscapeRadius;
 AR = PixelWidth*5.0; /*   radius of circle around attractor ZA = target set for attracting points */
 AR2 = AR*AR;
 ARn = GiveNormalizedParameters(); // compute normalized AR = ARn 
 
 //iMaxN = 10000*iMax;

 iXmax = iWidth; // range(iX) = [0, iXmax-1] 
 iYmax = iHeight; // range(iY) = [0, iYmax -1]




return 0; 
}

/* --------------------------------------------------------------------------------------------------------- */

int main(){

int i0;
int i1;

iPeriodChild =2;
setup();


i0 =  GiveAttractiveTime(ZCr); 
printf(" i0 = %d \n", i0);
i1 = GiveAttractiveTime( ZCr+I*PixelHeight*5.0); 
printf(" i1 = %d \n", i1);

  
 
 return 0;

}
 
