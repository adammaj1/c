
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

double m=6.65;  // = iXmax/iYmax

/* iXmax/iYmax = 1 */
unsigned int iSide = 1000; /* side of rectangle in pixels */
unsigned int iXmax; // ((int)(m*iSide)) /* height of image in pixels */
unsigned int iYmax ; //= iSide;
unsigned int iLength; // = (iXmax*iYmax) /* number of pixels */
/*  world ( double) coordinate */
double dKSide =  0.3333;
/*
--- a = periods 2-8 ---
KxMin = 0.0;
KxMax  = m*dKSide;
---- b = periods 8-14 ----
KxMin = m*dKSide + 0.03;
KxMax  = 2*m*dKSide + 0.06;
----- c = periods --
KxMin = 2*m*dKSide + 0.06;
KxMax = KxMin + m*dKSide;

*/
double KxMin ; // = 0.0;
double KxMax ; // =(m*dKSide);
double KyMin = 0.0;
double KyMax ; // = dKSide; 
/* (KxMax-KxMin)/(KyMax-KyMin)==iXmax/iYmax = 1 */

unsigned int IterationMax; // = (iXmax*100) /* proportional to resolution of picture */

double KPixelWidth; //=  ((KxMax-KxMin)/iXmax)
double KPixelHeight ;//= ((KyMax-KyMin)/iYmax)

double CDistanceMax ; //= PixelWidth; /* proportional to pixel size */


/* fc(z) = z*z + c */



double EscapeRadius = 2.0;  /* radius of circle around origin; its complement is a target set for escaping points */
double ER2 ; //= (EscapeRadius*EscapeRadius)


/* colors = shades of gray from 0=black  to 255=white */
unsigned char iExterior = 245; /* exterior of Julia set */
unsigned char iBoundary = 0; /* border , boundary*/
unsigned char iInterior = 180;







unsigned int f(unsigned int _iX, unsigned int _iY)
/* 
   gives position of point (iX,iY) in 1D array  ; uses also global variables 
   it does not check if index is good  so memory error is possible 
*/
{return (_iX + (iYmax-_iY-1)*iXmax );}





/*
http://mathr.co.uk/blog/2013-12-16_stretching_cusps.html
Claude Heiland-Allen

Figure 4.22 on pages 204-205 of The Science Of Fractal Images is presented with this description:

    A region along the cardioid is continuously blown up and stretched out, so that the respective segment of the cardioid becomes a line segment. ... Our blow-up factor is chosen accordingly to the result that all disks in Figure 4.22 have the same size.

transformation from line, thru half of the circle, to half of the cardioid 

There are 3 complex planes :
* k-plane ( where are lines) in Claude Heiland-Allen notation
* w-plane ( where are circles)
* c-plane ( for cardioid )

2 steps :
* from z plane go to w plane using w=fi(k)
* from w plane go to the  c plane using  c = gi(w)

so c = gi( fi ( k)) 

 <pre>
z(%i1) k:x+y*%i;
(%o1) %i*y+x
(%i2) fi(k):=(-%i-k)/(%i-k);
(%o2) fi(k):=−%i−k/%i−k
(%i3) gi(w):=w/2-w*w/4;
(%o3) gi(w):=w/2−(w*w)/4
(%i4) gfi(k):=gi(fi(k));
(%o4) gfi(k):=gi(fi(k))
(%i5) c:gfi(k)$
</pre>
How to compute c from k  without using CAS :  
<pre>
(%i6)  ratsimp(realpart(c));
(%o6) (y^4−4*y^3+(2*x^2+2)*y^2+(4−4*x^2)*y+x^4+6*x^2−3)/(4*y^4−16*y^3+(8*x^2+24)*y^2+(−16*x^2−16)*y+4*x^4+8*x^2+4)
(%i7) ratsimp(imagpart(c));
(%o7) −(2*x*y−2*x)/(y^4−4*y^3+(2*x^2+6)*y^2+(−4*x^2−4)*y+x^4+2*x^2+1)
</pre>

Check the gfi function with some known values :
<pre>
(%i8) gfi(0);
(%o8) −3/4
(%i12) gfi(%i/3);
(%o12) −2
</pre>


*/

double GiveCx(double x,double y)
{
  return  (y*y*y*y-4*y*y*y+(2*x*x+2)*y*y+(4-4*x*x)*y+x*x*x*x+6*x*x-3)
         /(4*y*y*y*y-16*y*y*y+(8*x*x+24)*y*y+(-16*x*x-16)*y+4*x*x*x*x+8*x*x+4) ;
}

double GiveCy(double x,double y)
{
  return -(2*x*y-2*x)/(y*y*y*y-4*y*y*y+(2*x*x+6)*y*y+(-4*x*x-4)*y+x*x*x*x+2*x*x+1);
}

//give distance between 2 pixels on c-plane
double GivePixelSize( int iX)
{
  double Cx1, Cy1, Cx0, Cy0;
  double Kx1, Ky1, Kx0, Ky0;
  double dx, dy;
  //double p;

  //p= 0.1 + iX/(iXmax);
  Ky0 = KyMin ; /*  */
  Ky1 = KyMin+KPixelHeight;
  Kx0 = KxMin + iX*KPixelWidth; // iX=0
  Kx1 = Kx0;
  //
  Cx0 = GiveCx(Kx0,Ky0);
  Cy0 = GiveCy(Kx0,Ky0);
  Cx1 = GiveCx(Kx1,Ky1);
  Cy1 = GiveCy(Kx1,Ky1);
  //
  dx = Cx1 - Cx0;
  dy = Cy1 - Cy0;
  //
  return sqrt(dx*dx+dy*dy);

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

 */

// boolean Escape time and DEM/M in one loop  
unsigned char GiveColor( double Kx, double Ky, int iMax, double CDistanceMax)
{ 

  double Cx, Cy; // C = Cx + Cy* I = point of parameter c-plane

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
  unsigned char color = iInterior;

  // mapping from k to c plane 
  // http://mathr.co.uk/blog/2013-12-16_stretching_cusps.html

  Cx = GiveCx(Kx,Ky);
  Cy = GiveCy(Kx,Ky); 

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
      if (absZ2  > ER2 ) { color=  iExterior; break ; } // exterior when escapes
      if (absdZ2 > 1e60) { return iInterior; } //  interior when derivative explodes 
      // z = fc(z) = z*z + c 
      Zy=2*Zx*Zy + Cy;
      Zx=Zx2-Zy2 +Cx;
      
     /* first derivative   zp = 2*z*zp  = xp + yp*i; */
      temp = 2*(Zx*dZx - Zy*dZy) + 1.0 ; /*  */ 
      dZy = 2*(Zx*dZy + Zy*dZx); 
      dZx = temp;

      // abs 
      Zx2 = Zx*Zx;      
      Zy2 = Zy*Zy;
      absZ2= Zx2 + Zy2; // nz = x*x + y*y;
      absdZ2= dZx*dZx + dZy*dZy; // 
    };

  if (i<iMax) 
   { 
     temp = sqrt(absZ2);
     distance = 2.0*temp*log(temp)/sqrt(absdZ2); //
     if (distance > 0.0) 
          { 
            if (distance < (CDistanceMax/1000.0) ) color = iBoundary ; // bondary = black
            // near boundary = shades of gray
            else if (distance < CDistanceMax/100 ) color = (iBoundary + (int)(255*distance/CDistanceMax)) % 255; 
        }    
    }
   // if (nz < bailout) return 1; //still not escaping after iteration , rare
  // if (absdZ2 < absZ2) color= iExterior;  //includes escaping through 0 // som eiterior points are coded as exterior = error
    
  return color; 
} //0.000007
 


/* --------------------------------------------------------------------------------------------------------- */

int main(){

  unsigned int iX,iY, /* indices of 2D virtual array (image) = integer coordinate */
    i; /* index of 1D array  */
  double Kx, Ky; // K is a coordinate from K-plane
  double CDistanceMax;  //  Max Distance on c-plane = parameter plane with Mandelbrot set 
  //int  ExtLastIteration; //
  unsigned char color;


  printf(" Setup \n");
  iXmax= (int)(m*iSide); /* height of image in pixels */
  iYmax  = iSide;
  iLength = iXmax*iYmax; /* number of pixels */
  KxMin = 2*m*dKSide + 0.06;
  KxMax = KxMin + m*dKSide  + 0.03;
  KyMax = dKSide; 
  IterationMax = iXmax*10 ;
  KPixelWidth=  (KxMax-KxMin)/iXmax;
  KPixelHeight = (KyMax-KyMin)/iYmax;
  ER2 = EscapeRadius*EscapeRadius;

  /* dynamic 1D array for colors ( shades of gray ) */
  unsigned char *data;
  data = malloc( iLength * sizeof(unsigned char) );
  if (data == NULL )
    {
      fprintf(stderr," Could not allocate memory");
      return 1;
    }
    
  printf(" compute color \n");
  for(iY=0;iY<iYmax;++iY){ 
    Ky=KyMin + iY*KPixelHeight; /*  */
    printf("row %u from %u \n",iY, iYmax);    
    for(iX=0;iX<iXmax;++iX){ 
      Kx=KxMin + iX*KPixelWidth;
      CDistanceMax = GivePixelSize( iX); 
      color=GiveColor(Kx, Ky, IterationMax, CDistanceMax);
      
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
      data[i]=color;  /* change the color */
	     
    }
  }


  

   

  printf(" save  data array to the pgm file \n");
  unsigned int length = 30;
  char filename[length] ;
  snprintf(filename, length, "%fe%u%s", EscapeRadius,iXmax, ".pgm");
  char *comment="#  ";/* comment should start with # */
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  /* save image to pgm file  */ 
  FILE * fp;     
  fp = fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P5\n %s\n %u\n %u\n %u\n",comment,iXmax,iYmax,MaxColorComponentValue);  /*write header to the file*/
  fwrite(data,iLength,1,fp);  /*write image data bytes to the file in one step */
  printf("File %s saved. \n", filename);
  fclose(fp);

  /* ---------- text file  -------------------------------------*/
  char tfilename[length] ;
  printf(" save  graphic data to the text file \n");
  snprintf(tfilename, length, "%fe%u%s",EscapeRadius, iXmax, ".txt");
   /* save info  to the file  */      
  fp = fopen(tfilename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"IterationMax = %d \n", IterationMax);
  fprintf(fp,"EscapeRadius = %f \n", EscapeRadius);
  fprintf(fp,"CDistanceMax = %.15f \n",  CDistanceMax);

  fprintf(fp,"\n" );  /* */
  fprintf(fp,"K plane : \n" );  /* */
  fprintf(fp,"dKSide = %f \n", dKSide);
  fprintf(fp,"\n" );  /* */
  //
  fprintf(fp,"Corners of rectangle on K-plane : \n ");  /* */
  fprintf(fp," KxMin = %f \n", KxMin);  /* */
  fprintf(fp," KxMax = %f \n", KxMax);  /* */
  fprintf(fp," KyMin = %f \n", KyMin);  /* */
  fprintf(fp," KyMax = %f \n", KyMax);  /* */
  fprintf(fp,"\n" );  /* */
  //
  fprintf(fp,"corners of strip on C-plane: \n ");  /* */
  fprintf(fp," CxMin = %f \n", GiveCx(KxMin,KyMin));  /* */
  fprintf(fp," CxMax = %f \n", GiveCx(KxMax,KyMin));  /* */
  fprintf(fp," CyMin = %f \n", GiveCx(KxMin,KyMax));  /* */
  fprintf(fp," CyMax = %f \n", GiveCx(KxMax,KyMax));  /* */
  fprintf(fp,"\n" );
  //
  printf("File %s saved. \n", tfilename);
  fclose(fp);


  printf(" allways free memory \n ");
  free(data);
  
  

  return 0;
}

