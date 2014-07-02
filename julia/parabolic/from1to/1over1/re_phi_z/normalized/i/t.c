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
 
 
  gcc t.c -lm -Wall -o2
  gcc t.c -lm -Wall -march=native
  time ./a.out
 
 
 convert -resize 1000x1000  10000001.002100000.pgm 002100000.png
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

/* iXmax/iYmax = 1 */
unsigned int iXmax = 2000; /* height of image in pixels */
unsigned int iYmax = 2000;

unsigned int iLength; 
/* fc(z) = z*z + c */
unsigned int  denominator= 1; /* denominator of internal angle */
//double CxMin = 0.000;
double Cx =  0.25; /* C = Cx + Cy*i */
double Cy =  0.0;

int IterationMax;
//#define alfa (1-sqrt(1-4*Cx))/2 /* attracting or parabolic fixed point z = alfa */
double AR ;  /* PixelWidth*1.5   radius of circle around attractor ZA = target set for attracting points */
double AR2; //= (AR*AR);


/* color */
//unsigned char color[]={255,231,123,99}; /* shades of gray used in image */
const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
unsigned char iExterior = 245;
unsigned char iInteriorOdd = 200;
unsigned char iInteriorEven = 255;
unsigned char iPetal = 100; // ???


 
 
  unsigned int iX,iY, /* indices of 2D virtual array (image) = integer coordinate */
    i; /* index of 1D array  */
   
  const double dSide = 1.5;
 double PixelWidth; //=(ZxMax-ZxMin)/iXmax;
  double PixelHeight; //=(ZyMax-ZyMin)/iYmax; 







unsigned char GiveColor(double Zx0, double Zy0,double Cx, double Cy, int iMax, double ER2, double AR2, double ZAx, double ZAy)
 {
int i;
  double Zx, Zy;
  double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
  double d2, dX, dY; /* distance from z to Alpha  */

  Zx = Zx0; /* initial value of orbit  */
  Zy = Zy0;
  Zx2 = Zx*Zx;
  Zy2 = Zy*Zy;
  dX = Zx-ZAx;
  dY = Zy-ZAy;
  d2 = dX*dX+dY*dY; // d2= d*d where d is a distance to fixed point alfa ZA=ZAx+ZAy*i


  for (i=0; i<iMax; i++)
    { 
      if (Zx2+Zy2>ER2) return iExterior; // escaping to infinity = exterior
      // attracting to alfa = falling into target set around alfa
      // color depends on ???
      //  
      if (d2< AR2) { if (i % 2) return iInteriorOdd; 
                             else return iInteriorEven;} 
      //  forward iteration of z under fc(z) = z^2 + c
      Zy=2*Zx*Zy + Cy;
      Zx=Zx2-Zy2 +Cx;
      //
      Zx2=Zx*Zx;
      Zy2=Zy*Zy;
      //
      dX = Zx-ZAx;
      dY = Zy-ZAy;
      d2 = dX*dX+dY*dY;
    };
  return iPetal; // ? never happens ??? 

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
 increse AR s

 here z_cr = 0.0
 changes iterationMax and AR 
*/
int GiveNormalizedParameters(double Cx, double Cy, double alfax)
{

 int i=0;
// initiall point of iteration is a critical point of f 
double Zx=0.0;
double  Zy=0.0;
  double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
 double  d; /* distance from z to Alpha  */


do
    { 
     //  forward iteration of z under fc(z) = z^2 + c
      Zy=2*Zx*Zy + Cy;
      Zx=Zx2-Zy2 +Cx;
      //
      Zx2=Zx*Zx;
      Zy2=Zy*Zy;
       d = alfax-Zx;
      i+=1; 
          }
  while (d>PixelWidth);
  // first output using global variable 
  AR = d;
  AR2= AR*AR;
 // second output
 return i*100; // not enough for other interior points then critical point 

}





 
/* gives position of point (iX,iY) in 1D array  ; uses also global variables */
unsigned int f(unsigned int _iX, unsigned int _iY)
{return (_iX + (iYmax-_iY-1)*iXmax );}


// save data array to pgm file 
int SavePGMFile(double Cx, unsigned char data[])
{
FILE * fp;
  char name [10]; /* name of file */
  sprintf(name,"%10.9f", Cx); /*  */
  char *filename =strcat(name,".pgm");
  char *comment="# ";/* comment should start with # */
  /* save image to the pgm file  */      
  fp= fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P5\n %s\n %u %u\n %u\n",comment,iXmax,iYmax,MaxColorComponentValue);  /*write header to the file*/
  fwrite(data,iLength,1,fp);  /*write image data bytes to the file in one step */
  printf("File %s saved. \n", filename);
  fclose(fp);
  return 0;
}



 
/* --------------------------------------------------------------------------------------------------------- */
 
int main(){
  /* world ( double) coordinate = dynamic plane = z-plane */
  double ZxMin=-dSide;
  double ZxMax=dSide;
  double ZyMin=-dSide;
  double ZyMax=dSide;
   
  /* */
  double Zx, Zy;    /* Z=Zx+Zy*i   */
  double alfax; // define alfa (1-sqrt(1-4*Cx))/2 /* attracting or parabolic fixed point z = alfa */
  //double complex ZA;  /* atractor ZA = ZAx + ZAy*i */
  /* */
 
  const double EscapeRadius=2.0; /* radius of circle around origin; its complement is a target set for escaping points */
  double ER2; //=EscapeRadius*EscapeRadius;
 
  
  
  
  /* sobel filter is used to finfd boundaries of level sets */
  unsigned char G, Gh, Gv; 

 PixelWidth=(ZxMax-ZxMin)/iXmax;
   PixelHeight=(ZyMax-ZyMin)/iYmax;
  ER2=EscapeRadius*EscapeRadius;
 
  iLength = iXmax*iYmax;/* length of array in bytes = number of bytes = number of pixels of image * number of bytes of color */
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
 
 
 alfax = (1-sqrt(1-4*Cx))/2 ; /* parabolic fixed point z = alfa */
 IterationMax = GiveNormalizedParameters(Cx, Cy, alfax); // iterationMax and AR
 
  
  
   
   
 
 
 
 
 printf(" fill the data array \n");
  for(iY=0;iY<iYmax;++iY){ 
    Zy=ZyMin + iY*PixelHeight; 
    if (fabs(Zy)<PixelHeight/2) Zy=0.0; 
   printf(" row %u from %u \n",iY, iYmax);    /* info */
    for(iX=0;iX<iXmax;++iX){ 
      Zx=ZxMin + iX*PixelWidth;
      
      i= f(iX,iY);  /* compute index of 1D array from indices of 2D array */
      data[i]= GiveColor(Zx, Zy, Cx, Cy, IterationMax, ER2, AR2, alfax, 0.0);

        
        //if (Zx>=0 && Zx <= 0.5 && (Zy > 0 ? Zy : -Zy) <= 0.5 - Zx) data[i]=255-data[i]; // show petal

         
      /*  if (Zx>0 && Zy>0) data[i]=255-data[i];    check the orientation of Z-plane by marking first quadrant */




    }
  }
 
 
 printf(" find boundaries of level sets in data array using  Sobel filter\n");   
 
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
  printf(" save  arrays to the pgm files \n");
  SavePGMFile( IterationMax+0.0+AR, data); // level sets 
  SavePGMFile( IterationMax+1.0+AR, edge); // Re(phi(z)) is integer number ( Z )

  
 
  /* --------------free memory ---------------------*/
  free(data);
  free(edge);
  /* ------------ info -------------------*/
  printf("PixelWidth= %.10f \n",PixelWidth);
  printf(" AR = %.10f \n", AR);
  printf(" alfax = %f \n", alfax);
 
 
  return 0;
}

