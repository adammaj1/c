
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
  gcc m.c -std=c99 -lm -Wall -march=native
  to run ( Linux console) :
  ./a.out


 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double m=6.65;  // = iXmax/iYmax

/* iXmax/iYmax = 1 */
unsigned int iSide = 500; /* side of rectangle in pixels */
unsigned int iXmax; // ((int)(m*iSide)) /* height of image in pixels */
unsigned int iYmax ; //= iSide;
unsigned int iLength; // = (iXmax*iYmax) /* number of pixels */
/*  world ( double) coordinate */
double dKSide =  0.3333;
double KxMin ;//= 0.0;
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
unsigned char iInterior = 230;



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
int GivePeriod(double Cx,double Cy, int Iteration_Max, double precision)
{  
  
 
  double Zx2, Zy2, /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
         ZPrevieousX,ZPrevieousY,
         ZNextX,ZNextY;
  
     int Iteration,I;
     int result=0; 
     unsigned int k=0; // index of the array 


   /*  */
  double *orbit;

   orbit = malloc(2*Iteration_Max*sizeof(double));
  if (orbit==NULL)
    {
      fprintf(stderr," Could not allocate memory");
      return 1;
    };
  
 
   //double orbit[Iteration_Max+1][2]; /* array elements are numbered from 0 to length-1 */    
 
  /* starting point is critical point  */
   ZPrevieousX=0.0;
   ZPrevieousY=0.0;
   orbit[2*k]=0.0;
   orbit[2*k+1]=0.0;  
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
            orbit[2*Iteration]=ZNextX;
            orbit[2*Iteration+1]=ZNextY;   
            
        }; 
    /* here iteration=IterationMax+1 but last element of orbit has number IterationMax */    
     for(I=2*Iteration_Max-1;I>0;I--) 
      if (SameComplexValue(orbit[2*Iteration_Max-1],orbit[2*Iteration_Max],orbit[2*I],orbit[2*I+1],precision))
        result=Iteration_Max-I;

  free(orbit);
  return result;
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

//give distance between 2 pixels on c-plane
double GiveCDistanceMax( int iX)
{
  double Cx1, Cy1, Cx0, Cy0;
  double Kx1, Ky1, Kx0, Ky0;
  double dx, dy;
  double p;

  p= 0.5 + iX/(iXmax);
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
  return sqrt(dx*dx+dy*dy)*p;

}

/* --------------------------------------------------------------------------------------------------------- */

int main(){

  unsigned int iX,iY, /* indices of 2D virtual array (image) = integer coordinate */
    i; /* index of 1D array  */
  double Kx, Ky; // K is a coordinate from K-plane
  double Cx,Cy;  //  c is a coordinate from c-plane = parameter plane with Mandelbrot set 
  int  ExtLastIteration; //     period;


  
  iXmax= (int)(m*iSide); /* height of image in pixels */
  iYmax  = iSide;
  iLength = iXmax*iYmax; /* number of pixels */
  KxMin = m*dKSide - 0.05;
  KxMax  = 2*KxMin;
  KyMax = dKSide; 
  IterationMax = iXmax*10 ;
  KPixelWidth=  (KxMax-KxMin)/iXmax;
  KPixelHeight = (KyMax-KyMin)/iYmax;
  ER2 = EscapeRadius*EscapeRadius;

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
 
  
  printf(" find components of Mandelbrot set and save them to the data array \n");
  for(iY=0;iY<iYmax;++iY){ 
    Ky=KyMin + iY*KPixelHeight; /*  */
    if (fabs(Cy)<KPixelHeight/2) Cy=0.0; /* use it for interior , not for boundary  */
    printf("Period row %u from %u \n",iY, iYmax);    
    for(iX=0;iX<iXmax;++iX){ 
      Kx=KxMin + iX*KPixelWidth;
      Cx = GiveCx(Kx,Ky);
      Cy = GiveCy(Kx,Ky); 
      ExtLastIteration = GiveExtLastIteration(Cx, Cy, IterationMax, ER2);
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
      if (  ExtLastIteration <IterationMax ) 
	 data[i]=iExterior;
         else  data[i]=iInterior;  /* interior */
	  
	
      /* if (Cx>0 && Cy>0) data[i]=255-data[i];    check the orientation of Z-plane by marking first quadrant */
    }
  }


  printf(" find boundaries of components of Mandelbrot set in data array using Sobel filter and save to edge array \n"); 
  unsigned char G, Gh, Gv;  /* sobel filter */ 
  for(iY=1;iY<iYmax-1;++iY){ 
    for(iX=1;iX<iXmax-1;++iX){ 
      Gv= data[f(iX-1,iY+1)] + 2*data[f(iX,iY+1)] + data[f(iX-1,iY+1)] - data[f(iX-1,iY-1)] - 2*data[f(iX-1,iY)] - data[f(iX+1,iY-1)];
      Gh= data[f(iX+1,iY+1)] + 2*data[f(iX+1,iY)] + data[f(iX-1,iY-1)] - data[f(iX+1,iY-1)] - 2*data[f(iX-1,iY)] - data[f(iX-1,iY-1)];
      G = sqrt(Gh*Gh + Gv*Gv);
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
      if (G==0) 
	{edge[i]=255;} /* background */
      else {edge[i]=iBoundary;}  /* boundary */
    }
  }

  printf(" find boundary of Mandelbrot set using DEM/M \n");
  for(iY=0;iY<iYmax;++iY){ 
    printf(" DEM row %u from %u \n",iY, iYmax); 
      
    for(iX=0;iX<iXmax;++iX){ 
      
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
      if ( data[i]==iExterior ) 
	{ CDistanceMax=GiveCDistanceMax( iX);
          Ky=KyMin + iY*KPixelHeight;
          //if (fabs(Cy)<KPixelHeight/2) Cy=0.0; /* use it for interior , not for boundary  */
	  Kx=KxMin + iX*KPixelWidth;
          Cx = GiveCx(Kx,Ky);
          Cy = GiveCy(Kx,Ky);
	  if (mDist(Cx,Cy,IterationMax+iY)<CDistanceMax) data[i]=iBoundary;}
	else data[i]= iInterior; 
      /* if (Cx>0 && Cy>0) data[i]=255-data[i];    check the orientation of Z-plane by marking first quadrant */
    }
  }

 

  printf(" copy components boundaries from edge to data array \n");
  for(iY=1;iY<iYmax-1;++iY){ 
    for(iX=1;iX<iXmax-1;++iX)
      {i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
	if (edge[i]==iBoundary) data[i]=iBoundary;}}
  
   

  /* ---------- graphic file  -------------------------------------*/
  printf(" save  data array to the pgm file \n");
  FILE * fp;
  char name [30]; /* name of file */
  sprintf(name,"d%fm%u",m,iXmax); /* result (is saved in i) but is not used */
  char *filename =strcat(name,".pgm");
  char *comment="#  ";/* comment should start with # */
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  /* save image to pgm file  */      
  fp = fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P5\n %s\n %u\n %u\n %u\n",comment,iXmax,iYmax,MaxColorComponentValue);  /*write header to the file*/
  fwrite(data,iLength,1,fp);  /*write image data bytes to the file in one step */
  printf("File %s saved. \n", filename);
  fclose(fp);

  /* ---------- text file  -------------------------------------*/
  printf(" save  data array to the pgm file \n");
  filename = strcat(name,".txt");
   /* save info  to the file  */      
  fp = fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"IterationMax = %d \n", IterationMax);
  fprintf(fp,"EscapeRadius = %f \n", EscapeRadius);

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
  printf("File %s saved. \n", filename);
  fclose(fp);


  /* --------------free memory ---------------------*/
  free(data);
  free(edge);
  
  

  return 0;
}

