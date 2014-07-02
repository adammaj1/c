
/*

  c console program, for CPU, one thread. numbers type : double 
  It can be compiled and run under Linux, windows, Mac 
  It needs gcc


draw :
* check 2 algorithms :
** binary escape time
** add boundary computed by DEM/M
* save it to the pgm file 

 
 
  



  -----------------------------------------
  1.pgm file code is  based on the code of Claudio Rocchini
  http://en.wikipedia.org/wiki/Image:Color_complex_plot.jpg
  create 8 bit color graphic file ,  portable gray map file = pgm 
  see http://en.wikipedia.org/wiki/Portable_pixmap
  to see the file use external application ( graphic viewer)
  I think that creating graphic can't be simpler
  ---------------------------
  2. first it creates data array which is used to store color values of pixels,
  fills tha array with data and after that writes the data from array to pgm file.
  It alows free ( non sequential) acces to "pixels"
    
  -------------------------------------------
  Adam Majewski   fraktal.republika.pl 
 
  

  to compile : 
  gcc c.c  -lm -Wall -march=native
  to run ( Linux console) :
  time ./a.out




convert h6.650000m6650.pgm -resize 800x120 c.png
convert b6.650000m6650.pgm -resize 1600x240 cbig.png
 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



/* iXmax/iYmax = 1 */
unsigned int iSide = 1000; /* side of rectangle in pixels */
unsigned int iXmax; // ((int)(m*iSide)) /* height of image in pixels */
unsigned int iYmax ; //= iSide;
unsigned int iLength; // = (iXmax*iYmax) /* number of pixels */
/*  world ( double) coordinate */
double dSide =  1.5;
double CxMin ; // = 0.0;
double CxMax ; // =(m*dSide);
double CyMin ; //= dSide;
double CyMax ; // = dSide; 
/* (CxMax-CxMin)/(CyMax-CyMin)==iXmax/iYmax = 1 */

unsigned int IterationMax; // = (iXmax*100) /* proportional to resolution of picture */

double PixelWidth; //=  ((CxMax-CxMin)/iXmax)
double PixelHeight ;//= ((CyMax-CyMin)/iYmax)

double CDistanceMax ; //= PixelWidth; /* proportional to pixel size */


/* fc(z) = z*z + c */



double EscapeRadius = 33.0;  /* radius of circle around origin; its complement is a target set for escaping points */
double ER2 ; //= (EscapeRadius*EscapeRadius)


/* colors = shades of gray from 0=black  to 255=white */
unsigned char iExterior = 255; /* exterior of Julia set */
unsigned char iBoundary = 0; /* border , boundary*/
unsigned char iInterior = 0;







unsigned int f(unsigned int _iX, unsigned int _iY)
/* 
   gives position of point (iX,iY) in 1D array  ; uses also global variables 
   it does not check if index is good  so memory error is possible 
*/
{return (_iX + (iYmax-_iY-1)*iXmax );}

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
int GivePeriod(double Cx,double Cy, int Iteration_Max, double precision, double Zp[2])
{  
 
 
  double Zx2, Zy2, /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
         ZPrevieousX,ZPrevieousY,
         ZNextX,ZNextY;
 
     int Iteration,
      I; 
     int  period = 0; // not periodic or period > Iteration_Max
  double orbit[Iteration_Max+1][2]; /* array elements are numbered from 0 to length-1 */    
 
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
 



//  color is proportional to ExtDistance between point c and nearest point in Mandelbrot set 
unsigned char GiveColor( double Cx, double Cy, int iMax, double DistanceMax)
{ 
  double ExtDistance ;
  unsigned char color ;
  int period;
  double Zp[2]; // periodic point 
  double IntDistance;
  
 
     ExtDistance = GiveDistance( Cx, Cy, iMax, DistanceMax);
     
     if (ExtDistance > 0.0) // exterior
        { if (ExtDistance<DistanceMax*333.0) // = clamp(ExtDistance, 0.0, 1.0) to remove level sets effect  !!!!!
             color = (int)(255.0*ExtDistance); // boundary and near boundary = shades of gray
             else color = iExterior;  // exterior far away from boundary 
        }
      else // interior  
       {
            period = GivePeriod(Cx, Cy, iMax, PixelWidth*10.0, Zp);
            //IntDistance = getInteriorLowerBound(period, Zp);
            // printf("period = %d ; IntDistance =  %f\n",period,  IntDistance);
            color = period*10 ; //255- (int)(255.0*IntDistance); 
            //color = iInterior;
       } 
    
  return color; 
} 
 

/* --------------------------------------------------------------------------------------------------------- */

int main(){

  unsigned int iX,iY, /* indices of 2D virtual array (image) = integer coordinate */
    i; /* index of 1D array  */
  double Cx, Cy; 
  
  //int  ExtLastIteration; //
  unsigned char color;


  printf(" Setup \n");
  iXmax = iSide; /* height of image in pixels */
  iYmax = iSide;
  iLength = iXmax*iYmax; /* number of pixels */
  CxMin = -2.25;
  CxMax = 0.75  ;
  CyMin = -dSide;
  CyMax = dSide; 
  IterationMax = 1000 ;
  PixelWidth=  (CxMax-CxMin)/iXmax;
  PixelHeight = (CyMax-CyMin)/iYmax;
  ER2 = EscapeRadius*EscapeRadius;

  /* dynamic 1D array for colors ( shades of gray ) */
  unsigned char *data;
  data = malloc( iLength * sizeof(unsigned char) );
  if (data == NULL )
    {
      fprintf(stderr," Could not allocate memory");
      return 1;
    }
    
  printf(" compute color for every pixel : scan c-plane \n");
  for(iY=0;iY<iYmax;++iY){ 
    Cy=CyMin + iY*PixelHeight; /*  */
    printf("row %u from %u \n",iY, iYmax);    
    for(iX=0;iX<iXmax;++iX){ 
      Cx=CxMin + iX*PixelWidth;
      color=GiveColor(Cx, Cy, IterationMax, PixelWidth);
      
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
      data[i]=color;  /* change the color */
	     
    }
  }


  

   

  printf(" save  data array to the pgm file \n");
  unsigned int length = 30;
  char filename[length] ;
  snprintf(filename, length, "%.0fe%u%s", EscapeRadius,iXmax  , ".pgm");
  char *comment="#  ";/* comment should start with # */
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  FILE * fp;     
  fp = fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P5\n %s\n %u\n %u\n %u\n",comment,iXmax,iYmax,MaxColorComponentValue);  /*write header to the file*/
  fwrite(data,iLength,1,fp);  /*write image data bytes to the file in one step */
  printf("File %s saved. \n", filename);
  fclose(fp);

  printf(" save  graphic data to the text file \n");
  char tfilename[length] ;
  snprintf(tfilename, length, "%.0fe%u%s",EscapeRadius, iXmax, ".txt");
  fp = fopen(tfilename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"IterationMax = %d \n", IterationMax);
  fprintf(fp,"EscapeRadius = %f ER2 = %f \n", EscapeRadius, ER2);
  fprintf(fp,"\n" );  /* */
  fprintf(fp,"C plane : \n" );  /* */
  fprintf(fp,"dWidth  = %f ; dHeight = %f \n",CxMax- CxMin, CyMax - CyMin );
  fprintf(fp,"PixelWidth = %f ; PixelHeight = %f \n", PixelHeight, PixelWidth);
  fprintf(fp," CxMin = %f \n", CxMin);  /* */
  fprintf(fp," CxMax = %f \n", CxMax);  /* */
  fprintf(fp," CyMin = %f \n", CyMin);  /* */
  fprintf(fp," CyMax = %f \n", CyMax);  /* */
  fprintf(fp," center of image : C = %f ; %f \n",CxMax -(CxMax-CxMin)/2.0, CyMax - (CyMax-CyMin)/2.0);  /* */
  fprintf(fp,"\n" );
  printf("File %s saved. \n", tfilename);
  fclose(fp);


  printf(" allways free memory \n ");
  free(data);
  
  

  return 0;
}

