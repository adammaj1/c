
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


  gcc i.c -lm -Wall -o2
  gcc i.c -lm -Wall -march=native
  ./a.out


 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>



const int iMax =2000;

/* iXmax/iYmax = 1 */
#define iXmax 2000 /* height of image in pixels */
#define iYmax 1000
/* fc(z) = z*z + c */

#define iPeriodChild 2 // Period of secondary component joined by root point with the parent component 
int iPeriodParent = 1; // main cardioid of Mandelbrot set

unsigned int denominator;/* denominator of internal angle */
  double InternalAngle;
double Cx;
double Cy; 
double complex c;
double complex ZA;  /* atractor ZA = ZAx + ZAy*i */
  double ZAx,ZAy;



#define AR 0.05  /* PixelWidth*1.5   radius of circle around attractor ZA = target set for attracting points */
#define AR2 AR*AR
double EscapeRadius=2.0; /* radius of circle around origin; its complement is a target set for escaping points */
double ER2; //=EscapeRadius*EscapeRadius;
//#define alfa (1-sqrt(1-4*Cx))/2 /* attracting or parabolic fixed point z = alfa */
//#define beta (1+sqrt(1-4*Cx))/2 /* repelling or parabolic fixed point z = beta */

/* color */
  unsigned char Colors[]={255,230,180}; /* shades of gray used in image */
  unsigned char color;
  unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  unsigned char iExterior = 245;
  unsigned char iInteriorL = 200;
  unsigned char iInteriorR = 220;



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



/* gives position of point (iX,iY) in 1D array  ; uses also global variables */
unsigned int f(unsigned int _iX, unsigned int _iY)
{return (_iX + (iYmax-_iY-1)*iXmax );}




unsigned char GiveColor(double Zx0, double Zy0)
{
  int i;
  double Zx, Zy;
  double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
  double d2, dX, dY; /* distance from z to Alpha  */


  Zx=Zx0; /* initial value of orbit  */
  Zy=Zy0;
  Zx2=Zx*Zx;
  Zy2=Zy*Zy;
  dX=Zx-ZAx;
  dY=Zy-ZAy;
  d2=dX*dX+dY*dY;
// Fall-in test : iterate z and check when point z is inside target set.
// Note that test is done after every f^{{period}} not after every f
  for (i=0; i<iMax ;i++)
    {
      if (Zx2+Zy2 > ER2) return iExterior;
      if (d2 < AR2) {if (i % 2) return iInteriorR; else return iInteriorL;} 
      Zy=2*Zx*Zy + Cy;
      Zx=Zx2-Zy2 +Cx;
      Zx2=Zx*Zx;
      Zy2=Zy*Zy;
      //
      dX=Zx-ZAx;
      dY=Zy-ZAy;
      d2=dX*dX+dY*dY;
    };

return 0;

}

 

/* --------------------------------------------------------------------------------------------------------- */

int main(){
  
 
    
  unsigned int iX,iY, /* indices of 2D virtual array (image) = integer coordinate */
    i, /* index of 1D array  */
    iLength = iXmax*iYmax;/* length of array in bytes = number of bytes = number of pixels of image * number of bytes of color */
  /* world ( double) coordinate = parameter plane*/
  const double ZxMin=-2.0;
  const double ZxMax=2.0;
  const double ZyMin=-1.0;
  const double ZyMax=1.0;
  double PixelWidth=(ZxMax-ZxMin)/iXmax;
  double PixelHeight=(ZyMax-ZyMin)/iYmax;
  /* */
  double Zx, Zy;    /* Z=Zx+Zy*i   */
  
  /* */
  
  
  
  
  int eLastIteration, iLastIteration;
  /* sobel filter */
  unsigned char G, Gh, Gv; 
  
  

  /* dynamic 1D arrays for colors ( shades of gray ) */
  unsigned char *data, *edge;
  data = malloc( iLength * sizeof(unsigned char) );
  edge = malloc( iLength * sizeof(unsigned char) );
  if (data == NULL || edge==NULL)
    {
      fprintf(stderr," Could not allocate memory");
      getchar(); 
      return 1;
    }
  else printf(" memory is OK\n");


   ER2=EscapeRadius*EscapeRadius;
 
  denominator = iPeriodChild;
  InternalAngle = 1.0/((double) denominator);
 
  c = GiveC(InternalAngle, 1.0, iPeriodParent) ; // internal radius= 1.0 gives root point = parabolic parameter   
 Cx=creal(c);
  Cy=cimag(c);
 
  
  ZA = GiveAlfaFixedPoint( c);
 ZAx=creal(ZA);
 ZAy = cimag(ZA);



  printf(" fill the data array \n");
  for(iY=0;iY<iYmax;++iY){ 
    Zy=ZyMin + iY*PixelHeight; /*  */
    if (fabs(Zy)<PixelHeight/2) Zy=0.0; /*  */
    printf(" row %u from %u \n",iY, iYmax);    
    for(iX=0;iX<iXmax;++iX){ 
      Zx=ZxMin + iX*PixelWidth;
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
      color = GiveColor(Zx, Zy);
      data[i]= color; /* exterior */
      /*  if (Zx>0 && Zy>0) data[i]=255-data[i];    check the orientation of Z-plane by marking first quadrant */
    }
  }


  printf(" find boundaries in data array using  Sobel filter\n");   

  for(iY=1;iY<iYmax-1;++iY){ 
    for(iX=1;iX<iXmax-1;++iX){ 
      Gv= data[f(iX-1,iY+1)] + 2*data[f(iX,iY+1)] + data[f(iX-1,iY+1)] - data[f(iX-1,iY-1)] - 2*data[f(iX-1,iY)] - data[f(iX+1,iY-1)];
      Gh= data[f(iX+1,iY+1)] + 2*data[f(iX+1,iY)] + data[f(iX-1,iY-1)] - data[f(iX+1,iY-1)] - 2*data[f(iX-1,iY)] - data[f(iX-1,iY-1)];
      G = sqrt(Gh*Gh + Gv*Gv);
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
      if (G==0) {edge[i]=255;} /* background */
      else {edge[i]=0;}  /* boundary */
    }
  }

    printf(" copy boundaries from edge to data array \n");
    for(iY=1;iY<iYmax-1;++iY){ 
     for(iX=1;iX<iXmax-1;++iX)
      {i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
  	if (edge[i]==0) data[i]=0;}}


  /* ---------- file  -------------------------------------*/
  printf(" save  data array to the file \n");
  FILE * fp;
  char name [100]; /* name of file */
  i = sprintf(name,"B%2.9f",AR); /* result (is saved in i) but is not used */
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
  

 printf("parameter c = ( %f ; %f ) \n", Cx, Cy); 
 printf("alfa fixed point z = ( %f ; %f )  \n", creal(ZA), cimag(ZA)); 
  

  return 0;
}
