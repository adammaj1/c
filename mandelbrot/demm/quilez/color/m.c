
/*

  c console program, for CPU, one thread. numbers type : double 
  It can be compiled and run under Linux, windows, Mac 
  It needs gcc


draw :
* check 2 algorithms :
** binary escape time
** add boundary computed by DEM/M
* save it to the ppm file 

 
 
  



  -----------------------------------------
  1.pgm file code is  based on the code of Claudio Rocchini
  http://en.wikipedia.org/wiki/Image:Color_complex_plot.jpg
  create 24 bit color graphic file ,  portable pixel map file = ppm 
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
  gcc m.c  -lm -Wall -march=native
  to run ( Linux console) :
  time ./a.out




convert h6.650000m6650.pgm -resize 800x120 c.png
convert b6.650000m6650.pgm -resize 1600x240 cbig.png
 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


double g=100.0;
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
{return ((_iX + (iYmax-_iY-1)*iXmax)*3 );}

// ((iYmax-iY-1)*iXmax+iX)*3; 

/* based on Delphi function by Witold J.Janik */
void GiveRainbowColor(double position,unsigned char Color32[])
{
  /* if position > 1 then we have repetition of colors it maybe useful    */
 
  if (position>1.0){if (position-(int)position==0.0)position=1.0; else position=position-(int)position;}
 
 
 
 
  unsigned char nmax=6; /* number of color segments */
  double m=nmax* position;
 
  int n=(int)m; // integer of m
 
  double f=m-n;  // fraction of m
  unsigned char t=(int)(f*255);
 
 
 
switch( n){
   case 0: {
      Color32[0] = 255;
      Color32[1] = t;
      Color32[2] = 0;
       break;
    };
   case 1: {
      Color32[0] = 255 - t;
      Color32[1] = 255;
      Color32[2] = 0;
       break;
    };
   case 2: {
      Color32[0] = 0;
      Color32[1] = 255;
      Color32[2] = t;
       break;
    };
   case 3: {
      Color32[0] = 0;
      Color32[1] = 255 - t;
      Color32[2] = 255;
       break;
    };
   case 4: {
      Color32[0] = t;
      Color32[1] = 0;
      Color32[2] = 255;
       break;
    };
   case 5: {
      Color32[0] = 255;
      Color32[1] = 0;
      Color32[2] = 255 - t;
       break;
    };
    default: {
      Color32[0] = 255;
      Color32[1] = 0;
      Color32[2] = 0;
       break;
    };
 
}; // case
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
     distance = pow( g*distance, 0.25 );
                 
    }
   else distance=0.0; // interior 
   // if (nz < bailout) return 1; //still not escaping after iteration , rare
  // if (absdZ2 < absZ2) color= iExterior;  //includes escaping through 0 // som eiterior points are coded as exterior = error
    
  return distance; 
} //0.000007
 
// gives width in estimated distnce of near boundary strip araound Mandelbrot set
// true boundary is = DistanceMax*0.1
double GiveDistanceMax(double PixelWidth, int iSide)
{
  return PixelWidth*iSide*0.0333; // = 0.1
}


//  color is proportional to distance between point c and nearest point in Mandelbrot set 
int ComputeColor32( double Cx, double Cy, int iMax, double DistanceMax, unsigned char color32[3] )
{ 
  double distance ;
 
 
 
     distance = GiveDistance( Cx, Cy, iMax, DistanceMax);
     
     if (distance > 0.0)
        { 
             
              
         if (distance<DistanceMax) // boundary
             { // boundary
                     color32[0] = 255-(int)(255.0*distance); 
                     color32[1] = 255-(int)(255.0*distance); 
                     color32[2] = 255-(int)(255.0*distance); 
                 }
             
              
              else { // exterior far away from boundary
                     GiveRainbowColor(distance, color32); //  gradient of 32 rgb color 
                     
                  }
          
        }
      else { // interior
                     color32[0] = 0; 
                     color32[1] = 0; 
                     color32[2] = 0; 
                 }
    
  return 0; 
} 
 


/* --------------------------------------------------------------------------------------------------------- */

int main(){

  unsigned int iX,iY, /* indices of 2D virtual array (image) = integer coordinate */
    i; /* index of 1D array  */
  double Cx, Cy; 
  
  //int  ExtLastIteration; //
  unsigned char color32[3];


  printf(" Setup \n");
  iXmax = iSide; /* height of image in pixels */
  iYmax = iSide;
  iLength = iXmax*iYmax*3; /* number of pixels* number of bytes per color  */
  CxMin = -2.25;
  CxMax = 0.75  ;
  CyMin = -dSide;
  CyMax = dSide; 
  IterationMax = 1000 ;
  PixelWidth=  (CxMax-CxMin)/iXmax;
  PixelHeight = (CyMax-CyMin)/iYmax;
  ER2 = EscapeRadius*EscapeRadius;
  CDistanceMax = GiveDistanceMax( PixelWidth,  iSide);

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
      ComputeColor32(Cx, Cy, IterationMax, CDistanceMax, color32);
      
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
      /* change the color */
      data[i]   = color32[0];  
      data[i+1] = color32[1];
      data[i+2] = color32[2];
	     
    }
  }


  

   

  printf(" save  data array to the pgm file \n");
  unsigned int length = 30;
  char filename[length] ;
  snprintf(filename, length, "%fg%.0fe%u%s",g, EscapeRadius,iXmax  , ".ppm");
  char *comment="#  ";/* comment should start with # */
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  FILE * fp;     
  fp = fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P6\n %s\n %u\n %u\n %u\n",comment,iXmax,iYmax,MaxColorComponentValue);  /*write header to the file*/
  fwrite(data,iLength,1,fp);  /*write image data bytes to the file in one step */
  printf("File %s saved. \n", filename);
  fclose(fp);

  printf(" save  graphic data to the text file \n");
  char tfilename[length] ;
  snprintf(tfilename, length, "%fg%.0fe%u%s",g, EscapeRadius,iXmax  , ".txt");
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
  fprintf(fp,"g= %f used for  distance = pow( g*distance, 0.25 ); \n", g);
  printf("File %s saved. \n", tfilename);
  fclose(fp);


  printf(" allways free memory \n ");
  free(data);
  
  

  return 0;
}

