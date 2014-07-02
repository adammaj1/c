
/*

  c console program
  -----------------------------------------
  1.ppm file code is  based on the code of Claudio Rocchini
  http://en.wikipedia.org/wiki/Image:Color_complex_plot.jpg
  create 24 bit color graphic file ,  portable pixmap file = PPM 
  see http://en.wikipedia.org/wiki/Portable_pixmap
  to see the file use external application ( graphic viewer)
  I think that creating graphic can't be simpler
  ---------------------------
  2. first it creates data array which is used to store rgb color values of pixels,
  fills tha array with data and after that writes the data from array to pgm file.
  It alows free ( non sequential) acces to "pixels"
    
  -------------------------------------------
  Adam Majewski   fraktal.republika.pl 
 
  Sobel filter 
  Gh = sum of six values ( 3 values of matrix are equal to 0 ). Each value is = pixel_color * filter_coefficients 


gcc t.c -lm -Wall -O2
./a.out


 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>




/* iXmax/iYmax = 1 */
unsigned int iXmax= 1000; /* height of image in pixels */
unsigned int iYmax =1000;

 /* world ( double) coordinate = parameter plane*/
  const double ZxMin=-1.7;
  const double ZxMax=1.7;
  const double ZyMin=-1.7;
  const double ZyMax=1.7;
  double PixelWidth; //=(ZxMax-ZxMin)/iXmax;
  double PixelHeight; //=(ZyMax-ZyMin)/iYmax;




/* fc(z) = z*z + c */
 
double Cx= -0.75; /* C = Cx + Cy*i */
double Cy=  0.0;


  const int IterationMax=60;
    //IterationMaxBig= 1000001;
  int eLastIteration; //, iLastIteration;
  /* sobel filter */
  unsigned char G, Gh, Gv; 
  /* color */
  //unsigned char color[]={255,230,180}; /* shades of gray used in image */
  unsigned char iExterior = 255;
  unsigned char iInterior = 230;
  unsigned char iBoundary = 0;
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  

  /* dynamic 1D arrays for colors ( shades of gray ) */
  unsigned char *data, *edge;
 



















/* escape time to infinity */
int GiveExtLastIteration(double _Zx0, double _Zy0,double C_x, double C_y, int iMax, double _ER2)
{ 
  int i;
  double Zx, Zy;
  double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
  Zx=_Zx0; /* initial value of orbit  */
  Zy=_Zy0;
  Zx2=Zx*Zx;
  Zy2=Zy*Zy;
  for (i=0;i<iMax && ((Zx2+Zy2)<_ER2);i++)
    {
      Zy=2*Zx*Zy + C_y;
      Zx=Zx2-Zy2 +C_x;
      Zx2=Zx*Zx;
      Zy2=Zy*Zy;
    };
  return i;
}


int SetPixel(int iX, int iY , int color, unsigned char *array){


array[((iYmax-iY-1)*iXmax+iX)]=color;
return 0;
}

/*
 http://www.cs.unc.edu/~mcmillan/comp136/Lecture6/Lines.html 
Bresenham algorithm

*/

 int DrawLine(int x0, int y0, int x1, int y1, unsigned char color, unsigned char *array)
    {
        
        int dy = y1 - y0;
        int dx = x1 - x0;
        int stepx, stepy;

        if (dy < 0) { dy = -dy;  stepy = -1; } else { stepy = 1; }
        if (dx < 0) { dx = -dx;  stepx = -1; } else { stepx = 1; }
        dy <<= 1;                                                  // dy is now 2*dy
        dx <<= 1;                                                  // dx is now 2*dx

        SetPixel( x0, y0,color, array); /* first pixel */
        if (dx > dy) {
            int fraction = dy - (dx >> 1);                         // same as 2*dy - dx
            while (x0 != x1) {
                if (fraction >= 0) {
                    y0 += stepy;
                    fraction -= dx;                                // same as fraction -= 2*dx
                }
                x0 += stepx;
                fraction += dy;                                    // same as fraction -= 2*dy
                SetPixel(x0, y0,color,array);
            }
        } else {
            int fraction = dx - (dy >> 1);
            while (y0 != y1) {
                if (fraction >= 0) {
                    x0 += stepx;
                    fraction -= dy;
                }
                y0 += stepy;
                fraction += dx;
                SetPixel(x0, y0,color,array);
            }
        }
   return 0;
    }



int DrawLinesD(double complex z0, double complex z1, unsigned char color, unsigned char *array)
{
 unsigned int iX0, iY0; // pixel coordinate
   unsigned int iX1, iY1; // pixel coordinate
  double dX0,dY0;

  dX0=creal(z0);
  dY0=cimag(z0);
  
//  if z  is in window then draw a line from z0 to z1 = part of ray 
          if (ZxMin <= dX0 && dX0 < ZxMax && ZyMin <= dY0 && dY0 < ZyMax)
            {
              // convert to pixel coordinates = translate from world to screen coordinate 
              iX0=(dX0-ZxMin)/PixelWidth;
              iY0=(dY0-ZyMin)/PixelHeight; /*  */
              iX1=(creal(z1)-ZxMin)/PixelWidth;
              iY1=(cimag(z1)-ZyMin)/PixelHeight; /*  */
              DrawLine(iX0, iY0, iX1, iY1, color, array);
              // symmetrical point
              iX0=(-dX0-ZxMin)/PixelWidth;
              iY0=(-dY0-ZyMin)/PixelHeight; /*  */
              iX1=(-creal(z1)-ZxMin)/PixelWidth;
              iY1=(-cimag(z1)-ZyMin)/PixelHeight; /*  */
              DrawLine(iX0, iY0, iX1, iY1, color, array);
             }
           
return 0;
}


/* gives position of point (iX,iY) in 1D array  ; uses also global variables */
unsigned int f(unsigned int _iX, unsigned int _iY)
{return (_iX + (iYmax-_iY-1)*iXmax );}


/* 
 principal square  root of complex number 
 http://en.wikipedia.org/wiki/Square_root

 z1= I;
  z2 = root(z1);
  printf("zx  = %f \n", creal(z2));
  printf("zy  = %f \n", cimag(z2));
*/
double complex root(double complex z)
{ 
  double x = creal(z);
  double y = cimag(z);
  double u;
  double v;
  double r = sqrt(x*x + y*y); 
  
  v = sqrt(0.5*(r - x));
  if (y < 0) v = -v; 
  u = sqrt(0.5*(r + x));
 return u + v*I;
}

double complex preimage(double complex z1, double complex z2,  double complex c)
{ 
  double complex zPrev;

  zPrev = root(creal(z1) - creal(c) + ( cimag(z1) - cimag(c))*I);
  // choose one of 2 roots 
  if (creal(zPrev)*creal(z2) + cimag(zPrev)*cimag(z2) > 0) 
        return zPrev ;     //  u+v*i
       else return -zPrev; // -u-v*i
}


// -------------------------------------------------------------------
// based on backray function from program Mandel by Wolf Jung  
// http://www.mndynamics.com/indexp.html
//
// draws dynamic external ray using backward iteration
// This function only works for periodic or preperiodic angles.
// You must determine the period n and the preperiod k before calling this function.
double complex DrawRays(unsigned int numerator,
                        unsigned int denominator, // of external angle 
		       int n, //period of ray's angle under doubling map
		       int k, // preperiod
                       int iterMax,
                       double complex c,
                       unsigned char color, 
                       unsigned char *array
                          )
{
  const double R = 4; // very big radius = near infinity
  int j; // number of ray 
  int iter; // index of backward iteration
  
 
  
  double t; // external angle in turns
  double complex zNm1=0; // zNm1 = f^(N-1)(zN)





   /* dynamic 1D arrays for coordinates ( x, y) of points with the same R on preperiodic and periodic rays  */
  double *RayXs, *RayYs;
  int iLength = n+k+2; // length of arrays ?? why +2
  //  creates arrays :  RayXs and RayYs  and checks if it was done
  RayXs = malloc( iLength * sizeof(double) );
  RayYs = malloc( iLength * sizeof(double) );
  if (RayXs == NULL || RayYs==NULL)
    {
      fprintf(stderr,"Could not allocate memory");
      getchar(); 
      return 1;
    }
  
 t=(double) numerator/denominator; // angle of external ray

 //  starting points on preperiodic and periodic rays 
 //  with angles t, 2t, 4t...  and the same radius R
  for (j = 0; j < n + k; j++)
  { // z= R*exp(2*Pi*t)
    RayXs[j] = R*cos((2*M_PI)*t); 
    RayYs[j] = R*sin((2*M_PI)*t);
    t *= 2; // t = 2*t
    if (t > 1) t--; // t = t modulo 1 
  }
 // z[k] is n-periodic. So it can be defined here explicitly as well.
  RayXs[n+k] = RayXs[k]; 
  RayYs[n+k] = RayYs[k];


  //   backward iteration of starting points z on each ray
  for (iter = -10; iter <= iterMax; iter++)
    { 
     	
        for (j = 0; j < n+k; j++) // period +preperiod
         { // zNm1 = sqrt(z-c)   backward iteration of fc  
		zNm1 = root(RayXs[j+1] - creal(c)+(RayYs[j+1] - cimag(c))*I ); // 
                               
		// choose one of 2 roots: ZNm1 or -zNm1 
		if (creal(zNm1)*RayXs[j] + cimag(zNm1)*RayYs[j] <= 0) zNm1=-zNm1;
               

              // draw part of the ray for 
             DrawLinesD( RayXs[j]+ I*RayYs[j], zNm1, color, array);
              // save last point on ray
              RayXs[j] = creal(zNm1); 
              RayYs[j] = cimag(zNm1);    
         } // for j ...

          //RayYs[n+k] cannot be constructed as a preimage of RayYs[n+k+1]
          RayXs[n+k] = RayXs[k]; 
          RayYs[n+k] = RayYs[k];
          
     
  }

 
 
 // free memmory
 free(RayXs);
 free(RayYs);

return  creal(zNm1) + cimag(zNm1)*I; // return last point or ray for angle t 
}







/* --------------------------------------------------------------------------------------------------------- */

int main(){
  
 
    
  unsigned int iX,iY, /* indices of 2D virtual array (image) = integer coordinate */
    i, /* index of 1D array  */
    iLength = iXmax*iYmax;/* length of array in bytes = number of bytes = number of pixels of image * number of bytes of color */
 
  /* */
  double Zx, Zy;    /* Z=Zx+Zy*i   */
 // double complex ZA;  /* atractor ZA = ZAx + ZAy*i */
 double complex landing;
  /* */
  
  const double EscapeRadius=2.0; /* radius of circle around origin; its complement is a target set for escaping points */
  double ER2=EscapeRadius*EscapeRadius;
  
 data = malloc( iLength * sizeof(unsigned char) );
  edge = malloc( iLength * sizeof(unsigned char) );
  if (data == NULL || edge==NULL)
    {
      fprintf(stderr," Could not allocate memory");
      getchar(); 
      return 1;
    }
  else printf(" memory is OK\n");


   
 PixelWidth =(ZxMax-ZxMin)/iXmax;
  PixelHeight=(ZyMax-ZyMin)/iYmax;
  
 
   printf(" fill the data array \n");
  for(iY=0;iY<iYmax;++iY){ 
    Zy=ZyMin + iY*PixelHeight; /*  */
    if (fabs(Zy)<PixelHeight/2) Zy=0.0; /*  */
    printf(" row %u from %u \n",iY, iYmax); // info tekst   
    for(iX=0;iX<iXmax;++iX){ 
      Zx=ZxMin + iX*PixelWidth;
      eLastIteration = GiveExtLastIteration(Zx, Zy, Cx, Cy, IterationMax, ER2 );
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
      if ( IterationMax != eLastIteration ) 
	data[i]=iExterior; /* exterior */
      else data[i]=iInterior;/* interior */
	
      /*  if (Zx>0 && Zy>0) data[i]=255-data[i];    check the orientation of Z-plane by marking first quadrant */
    }
  }


printf(" use external rays to find narrow parts of exterior near parabolic fixed point and its preimages\n");
  
landing = DrawRays(1,3, 2, 0, 200000, Cx + Cy*I, iExterior, data );
 printf("landing =  (%f , %f ) \n",creal(landing), cimag(landing) );






   printf(" find boundaries in data array using  Sobel filter\n");   

  for(iY=1;iY<iYmax-1;++iY){ 
    for(iX=1;iX<iXmax-1;++iX){ 
      Gv= data[f(iX-1,iY+1)] + 2*data[f(iX,iY+1)] + data[f(iX-1,iY+1)] - data[f(iX-1,iY-1)] - 2*data[f(iX-1,iY)] - data[f(iX+1,iY-1)];
      Gh= data[f(iX+1,iY+1)] + 2*data[f(iX+1,iY)] + data[f(iX-1,iY-1)] - data[f(iX+1,iY-1)] - 2*data[f(iX-1,iY)] - data[f(iX-1,iY-1)];
      G = sqrt(Gh*Gh + Gv*Gv);
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
      if (G==0) {edge[i]=255;} /* background */
      else {edge[i]=iBoundary;}  /* boundary */
    }
  }

  printf(" copy boundaries from edge to data array \n");
  for(iY=1;iY<iYmax-1;++iY){ 
    for(iX=1;iX<iXmax-1;++iX)
      {i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
	if (edge[i]==iBoundary) data[i]=iBoundary;}}



  /* ---------- file  -------------------------------------*/
  printf(" save  data array to the file \n");
  FILE * fp;
  char name [20]; /* name of file */
  i = sprintf(name,"%12.10f",PixelWidth); /* result (is saved in i) but is not used */
  char *filename =strcat(name,".pgm");
  char *comment="# C=0.2";/* comment should start with # */
  /* save image to the pgm file  */      
  fp= fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P5\n %s\n %u\n %u\n %u\n",comment,iXmax,iYmax,MaxColorComponentValue);  /*write header to the file*/
  fwrite(data,iLength,1,fp);  /*write image data bytes to the file in one step */
  printf("File %s saved. \n", filename);
  fclose(fp);


  /* --------------free memory ---------------------*/
  free(data);
  free(edge);

  printf("PixelWidth =  %f \n", PixelWidth);
  
  

  return 0;
}
