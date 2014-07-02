
/*

  c console program, for CPU, one thread. numbers type : double 
  It can be compiled and run under Linux, windows, Mac 
  It needs gcc

draw :
- components using period 
- find boundaries of components using sobel filter
- add boundary computed by DEM/M
- save it to the pgm file 

 
  



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
  gcc m.c -lm -Wall -march=native
  to run ( Linux console) :
  ./a.out


 gcc -g m.c -lm
 gdb ./a.out


 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



#define EscapeRadius 2.0 /* radius of circle around origin; its complement is a target set for escaping points */
#define ER2 (EscapeRadius*EscapeRadius)


/* colors = shades of gray from 0=black  to 255=white */
#define iExterior 245 /* exterior of Julia set */
#define iBoundary 0 /* border , boundary*/
#define iInterior 230

int m=5;  // = iXmax/iYmax

/* iXmax/iYmax = 1 */
#define iSide 1000 /* side of rectangle in pixels */
#define iXmax (m*iSide) /* height of image in pixels */
#define iYmax iSide
#define iLength (iXmax*iYmax) /* number of pixels */

#define IterationMax (iXmax*10) /* proportional to resolution of picture */

/*  world ( double) coordinate */
double dSide = 0.3333;
double KxMin = 0.0;
double KxMax ; // = m*dSide;
double KyMin = 0.0;
double KyMax ; //= dSide;
/* (KxMax-KxMin)/(KyMax-KyMin)==iXmax/iYmax = 1 */



double PixelWidth ; // =   ((KxMax-KxMin)/iXmax)
double PixelHeight; // = ((KyMax-KyMin)/iYmax)

double distanceMax; // = (1.5*PixelWidth) /* proportional to pixel size */


/* fc(z) = z*z + c */







/* escape time to infinity of function fc(z) = z*z + c */
int GiveExtLastIteration(double C_x, double C_y, int iMax, double _ER2)
{ 
  int i; /* iteration */
  double Zx, Zy; /* Z = Zx + Zy*i */
  double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
  Zx=0.0; /* initial value of orbit  */
  Zy=0.0;
  Zx2=Zx*Zx;
  Zy2=Zy*Zy;
  for (i=0;i<iMax && ((Zx2+Zy2)<_ER2);i++)
    {
      Zy=2*Zx*Zy + C_y;
      Zx=Zx2-Zy2 +C_x;
      Zx2=Zx*Zx;
      Zy2=Zy*Zy;
    };
  return i; /* last iteration */
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
int GivePeriod(double Cx, double Cy, int Iteration_Max, double precision)
{  
  
  
  double Zx2, Zy2, /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
         ZPrevieousX,ZPrevieousY,
         ZNextX,ZNextY;
  int Iteration,Iter;  
  double orbit[Iteration_Max+1][2]; /* array elements are numbered from 0 to length-1 */    
  printf("Iteration = B\n");
  /* starting point is critical point  */
   ZPrevieousX=0.0;
   ZPrevieousY=0.0;
   orbit[0][0]=0.0;
   orbit[0][1]=0.0;  
   Zx2=ZPrevieousX*ZPrevieousX;
   Zy2=ZPrevieousY*ZPrevieousY;
   /* iterate and save points for analysis */
   for (Iteration=1;Iteration<Iteration_Max+1 ;Iteration++)
        {   printf("Iteration = %d \n", Iteration);
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
     for(Iter=Iteration_Max-1;Iter>0;Iter--) 
      if (SameComplexValue(orbit[Iteration_Max][0],orbit[Iteration_Max][1],orbit[Iter][0],orbit[Iter][Iter],precision))
        return(Iteration_Max-Iter);
  return 0;
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


 */
 double mDist(double Cx, double Cy ,  int iter_max)
 { 
 int i;
 double x = 0.0, /* Z = x+y*i */
         y = 0.0, 
         /* Zp = xp+yp*1 = 1  */
         xp = 1, 
         yp = 0, 
         /* temporary */
         nz,  
         nzp,
         /* a = abs(z) */
         a; 
 for (i = 1; i <= iter_max; i++)
  { /* first derivative   zp = 2*z*zp  = xp + yp*i; */
    nz = 2*(x*xp - y*yp) +1.0 ; /* ?? */ 
    yp = 2*(x*yp + y*xp); 
    xp = nz;
    /* z = z*z + c = x+y*i */
    nz = x*x - y*y + Cx; 
    y = 2*x*y + Cy; 
    x = nz; 
    /* */
    nz = x*x + y*y; 
    nzp = xp*xp + yp*yp;
    if (nzp > 1e60 || nz > 1e60) break;
  }
 a=sqrt(nz);
 /* distance = 2 * |Zn| * log|Zn| / |dZn| */
 return 2* a*log(a)/sqrt(nzp); 
 }



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
* z-plane ( where are lines) = k-plane in Claude Heiland-Allen notation
* w-plane ( where are circles)
* c-plane ( for cardioid )

2 steps :
* from z plane go to w plane using w=fi(z)
* from w plane go to the  c plane using  c = gi(w)

so c = gi( fi ( z)) 

 z:x+y*%i;
fi(z):=(-%i-z)/(%i-z);
gi(w):=w/2-w*w/4;
gfi(z):=gi(fi(z));
c:gfi(z)$


(%i7) ratsimp(realpart(c));
(%o7) (y^4−4*y^3+(2*x^2+2)*y^2+(4−4*x^2)*y+x^4+6*x^2−3)/(4*y^4−16*y^3+(8*x^2+24)*y^2+(−16*x^2−16)*y+4*x^4+8*x^2+4)
(%i8) ratsimp(imagpart(c));
(%o8) −(2*x*y−2*x)/                                     (y^4−4*y^3+(2*x^2+6)*y^2+(−4*x^2−4)*y+x^4+2*x^2+1)
 

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



/* --------------------------------------------------------------------------------------------------------- */

int main(){

  unsigned int iX,iY, /* indices of 2D virtual array (image) = integer coordinate */
    i; /* index of 1D array  */
  
  double Kx, Ky; // K is a coordinate from K-plane
  double Cx,Cy;  //  c is a coordinate from c-plane = parameter plane with Mandelbrot set 
  int //period, 
   ExtLastIteration;
  
  /* color */
  //unsigned char ColorList[]={255,230,180};
  
  /* two dynamic 1D arrays for colors ( shades of gray ) */
  unsigned char *data, *edge;
  data = malloc( iLength * sizeof(unsigned char) );
  edge = malloc( iLength * sizeof(unsigned char) );
  if (data == NULL || edge==NULL)
    {
      fprintf(stderr," Could not allocate memory");
      return 1;
    }
  else printf(" memory is OK\n");


 KxMax = m*dSide;
 printf("KxMax = %f \n",KxMax);
 KyMax = dSide;
 printf("KyMax = %f \n",KyMax);
 PixelWidth =   (KxMax-KxMin)/iXmax;
 printf("PixelWidth = %f \n",PixelWidth);
 PixelHeight = (KyMax-KyMin)/iYmax;
 printf("PixelHeight = %f \n",PixelHeight);
 distanceMax= 1.5*PixelWidth; /* proportional to pixel size */
 printf(" iLength = %u \n", iLength);
 printf("ER2 = %f \n", ER2);

         
 for(iY=0;iY<iYmax;++iY){ 
    Ky=KyMin + iY*PixelHeight; 
    //printf("Ky = %f \n",Ky);
    printf("Period row %u from %u \n",iY, iYmax);    
    for(iX=0;iX<iXmax;++iX){ 
         //printf("Period column %u from %u \n",iX, iXmax);    
         Kx=KxMin + iX*PixelWidth;
         //printf("Kx = %f \n",Kx);
         Cx = GiveCx(Kx,Ky);
         //printf("Cx = %f \n",Cx);
         Cy = GiveCy(Kx,Ky);  
         //printf("Cy = %0.15f \n",Cy);
         //printf(" IterationMax = %u \n", IterationMax);
         //printf("Iteration = A\n");
         //period = GivePeriod(Cx,Cy, IterationMax,  PixelWidth);
         ExtLastIteration = GiveExtLastIteration(Cx , Cy, IterationMax, ER2);
         //printf("Period = %u \n", period); 
         i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
         //printf("i = %u \n", i); 
         if ( ExtLastIteration <  IterationMax) 
	    data[i]=iExterior;
            else  data[i]=iInterior;  /* interior */
	  
     /* if (Cx>0 && Cy>0) data[i]=255-data[i];    check the orientation of Z-plane by marking first quadrant */
    }
  }






  printf(" find boundary of Mandelbrot set using DEM/M \n");
  for(iY=0;iY<iYmax;++iY){ 
    printf(" DEM row %u from %u \n",iY, iYmax); 
    Ky=KyMin + iY*PixelHeight; /*  */  
    for(iX=0;iX<iXmax;++iX){ 
      
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
      if ( data[i]==iExterior ) 
	{  
          Kx=KxMin + iX*PixelWidth;
          Cx = GiveCx(Kx,Ky);
          Cy = GiveCy(Kx,Ky);    
         if (mDist(Cx,Cy,IterationMax)<distanceMax) data[i]=iBoundary;}
	else data[i]= iInterior; 
      /* if (Cx>0 && Cy>0) data[i]=255-data[i];    check the orientation of Z-plane by marking first quadrant */
    }
  }

 

  printf(" copy components boundaries from edge to data array \n");
  for(iY=1;iY<iYmax-1;++iY){ 
    for(iX=1;iX<iXmax-1;++iX)
      {i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
	if (edge[i]==iBoundary) data[i]=iBoundary;}}
  
   

  /* ---------- file  -------------------------------------*/
  printf(" save  data array to the pgm file \n");
  FILE * fp;
  char name [10]; /* name of file */
  i = sprintf(name,"%um%u",m,iXmax); /* result (is saved in i) but is not used */
  char *filename =strcat(name,".pgm");
  char *comment="#  ";/* comment should start with # */
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  /* save image to the file  */      
  fp = fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P5\n %s\n %u\n %u\n %u\n",comment,iXmax,iYmax,MaxColorComponentValue);  /*write header to the file*/
  fwrite(data,iLength,1,fp);  /*write image data bytes to the file in one step */
  printf("File %s saved. \n", filename);
  fclose(fp);


  /* --------------free memory ---------------------*/
  free(data);
  free(edge);
  
  

  return 0;
}

