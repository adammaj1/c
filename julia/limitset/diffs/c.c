
 
/*
 
  c console program
  It can be compiled and run under Linux, windows, Mac 
  It needs gcc
 
  --------------------------------------
  draws critical orbit for f(z)=z*z+c 
  
 
 
 
  ------------------------------------------
  one can change :
  
  - n 
  - iSide ( width of image = iXmax = (iSide) 
  - NrOfCrOrbitPoints = ; // check rang of type for big numbers : iMax, i, NrOfCrOrbitPoints


 %lld and %llu for print long long int
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
  Here are 4 items : 
  1. complex plane Z with points z = zx+zy*i ( dynamic plane and world coordinate )
  2. virtual 2D array indexed by iX and iYmax ( integer or screen coordinate )
  3. memory 1D array data  indexed by i =f(iX,iY)
  4. pgm file 
 
 
 
  Adam Majewski   fraktal.republika.pl 
 
 
  to compile : 
  gcc c.c -lm -Wall -march=native
  to run ( Linux console) :
  time ./a.out



   p = 2 ; aver diff = 0.547566 
 p = 3 ; aver diff = 1.374846 
 p = 4 ; aver diff = 0.055426 
 p = 5 ; aver diff = 1.474369 
 p = 6 ; aver diff = 0.569561 
 p = 7 ; aver diff = 1.422716 
 p = 8 ; aver diff = 0.000688  !!!!!!
 p = 9 ; aver diff = 1.466179 
 p = 10 ; aver diff = 0.565243 
 p = 11 ; aver diff = 1.430001 
 p = 12 ; aver diff = 0.032035 
 p = 13 ; aver diff = 1.459237 
 p = 14 ; aver diff = 0.568260 
 p = 15 ; aver diff = 1.435870 
 p = 16 ; aver diff = 0.000406 
 p = 17 ; aver diff = 1.457347 
 p = 18 ; aver diff = 0.566840 
 p = 19 ; aver diff = 1.437839 





*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>





 
 
/* iXmax/iYmax =  */
const int iSide = 1000;
int iXmax ; /* width of image in pixels = (15*iSide); */
int iYmax ;
int iLength ;
/* */
const double ZxMin = -1.4;
const double ZxMax = 0.5;
const double ZyMin = -0.1;
const double ZyMax = 0.1;
/* (ZxMax-ZxMin)/(ZyMax-ZyMin)= iXmax/iYmax  */
 
 
double PixelWidth ;
double PixelHeight ;
 
 //unsigned int nMin = 0;
  //unsigned int nMax=7;
 unsigned int period=9;
 // unsigned int m;
 // check rang of type for big numbers : iMax, i, NrOfCrOrbitPoints
  unsigned long long int NrOfCrOrbitPoints ;
 
 
 
 
/* fc(z) = z*z + c */
/*   */
double Cx; /* C = Cx + Cy*i   */
double Cy ;
 
 
 
/* colors */
const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
const int iExterior = 245; /* exterior of Julia set */
const int iJulia = 0; /* border , boundary*/
const int iInterior = 230;
 
 
 
/* ----------------------- functions ---------------------------------------- */
 
 
unsigned int f(unsigned int iX, unsigned int iY)
/* 
   gives position of point (iX,iY) in 1D array  ; uses also global variables 
   it does not check if index is good  so memory error is possible 
*/
{return (iX + iY*iXmax );}
 
 
 
int DrawPoint( double Zx, double Zy, unsigned char data[])
{
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int index; /* index of 1D array  */ 
 
  
 
  iX = (int)((Zx-ZxMin)/PixelWidth);
  iY = (int)((Zy-ZyMin)/PixelHeight);
  index = f(iX,iY);
  
 
  data[index] = iJulia;  /* draw */
  return 0;
}
 

 /* save z to data file  */
double SaveZ2DataFile( double Zx, double Zy, unsigned int iMax)
{
  unsigned long long int i,j; /* nr of point of critical orbit */
  double Zx2,Zy2;
  double oldZx, oldZy;
  double diff;  
  double diffsum=0.0;
  double averagediff;
  char name [10]; /* name of file */
  sprintf(name,"data%u", iMax); /*  */
  char *filename =strcat(name,".txt");
  FILE * file;
  file = fopen(filename, "w");  
  if (file  == NULL)
    {
      fprintf(stderr, "Nie moge otworzyc %s\n",  filename);
      return 1;
    };
   
  fprintf( file,"# Cx = %10.10f Cy = %10.10f \n",Cx,Cy);
   /* titles  of columns */
    fprintf(file,"%s %s      %s \n","#","Zx","Zy");
 for (j=0; j<10; j++)
  {
   for (i=0; i<iMax; i++)
    {
      /* f(z) = z*z+c */
      Zy=2*Zx*Zy + Cy;
      Zx=Zx2-Zy2 + Cx;
 
      Zx2=Zx*Zx;
      Zy2=Zy*Zy;
      if (Zx2+Zy2>4) { printf("   point z = %f , %f escapes \n",Zx, Zy); break;}
      fprintf( file,"%10.10f %10.10f \t",Zx,Zy);
      
    }
    fprintf( file,"\n");
  }



 
  oldZx = Zx;
  oldZy = Zy;
  fprintf(file,"%s  %s \n","#","diff ");
  
  for (j=0; j<iMax; j++)
  {
   for (i=0; i<iMax; i++)
    {
      /* f(z) = z*z+c */
      Zy=2*Zx*Zy + Cy;
      Zx=Zx2-Zy2 + Cx;
      
      Zx2=Zx*Zx;
      Zy2=Zy*Zy;
      if (Zx2+Zy2>4) { printf("   point z = %f %f escapes \n",Zx, Zy); break;}
      
      
    }
    diff = sqrt((oldZx-Zx)*(oldZx - Zx) +(oldZy-Zy)*(oldZy-Zy));
    diffsum+=diff;
    fprintf( file,"%10.10f \t",diff);
    oldZx = Zx;
    oldZy = Zy;
   } 
  fprintf( file,"\n");
  averagediff  = diffsum/(1.0*iMax);
  fprintf( file,"# average diff = %10.10f \n",averagediff);  


  fclose( file);
  return averagediff;
}

 
// check rang of type for big numbers : iMax, i, NrOfCrOrbitPoints
double DrawLimitSet(unsigned long long int iMax,  unsigned char data[], unsigned int period )
{
  unsigned long long int i; /* nr of point of critical orbit */
  double Zx,Zy;
  double Zx2,Zy2; 
double averagediff;
 
  /* critical point z = 0 */
  Zx = 0.0;
  Zy = 0.0;
  DrawPoint(Zx,Zy,data);
  Zx2=Zx*Zx;
  Zy2=Zy*Zy;
  /* forward orbit of critical point  */
  for (i=1;i<=iMax ;i++)
    {
      /* f(z) = z*z+c */
      Zy=2*Zx*Zy + Cy;
      Zx=Zx2-Zy2 + Cx;
 
      Zx2=Zx*Zx;
      Zy2=Zy*Zy;
      if (Zx2+Zy2>4) { printf("   point z = %f , %f escapes \n",Zx, Zy); break;}
        printf(" %u ;  %llu \n",period, i); /* progres info */
    }

    for (i=1;i<=100 ;i++)
    {
      /* f(z) = z*z+c */
      Zy=2*Zx*Zy + Cy;
      Zx=Zx2-Zy2 + Cx;
 
      Zx2=Zx*Zx;
      Zy2=Zy*Zy;
      if (Zx2+Zy2>4) { printf("   point z = %f , %f escapes \n",Zx, Zy); break;}
      /* draws critical orbit */
      DrawPoint(Zx,Zy,data); 
      printf("  %llu \n",i); /* progres info */
    }
    /* save last z , maybe it will be used */
    averagediff=  SaveZ2DataFile(Zx,Zy, period);
  
  return averagediff;
 
}


int ClearArray( unsigned char data[] )
{
  unsigned int index; /* index of 1D array  */
  for(index=0;index<iLength-1;++index) 
                data[index]=iExterior;
  return 0;
}

int SaveArray2pgm(unsigned char data[], unsigned int n)
{

FILE * fp;
  char name [20]; /* name of file */
  sprintf(name,"data%u",n); /*  */
  char *filename =strcat(name,".pgm");
  char *comment="# C= ";/* comment should start with # */
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
  
  
 
  
 unsigned int p;
double diffs[2*period+1];
 
  /* unsigned int iX,iY;  indices of 2D virtual array (image) = integer coordinate */
  iXmax = iSide; /* height of image in pixels */
  iYmax = iSide;
  iLength = (iXmax*iYmax);
 
  PixelWidth =  ((ZxMax-ZxMin)/iXmax);
  PixelHeight = ((ZyMax-ZyMin)/iYmax);
 
 
  /* dynamic 1D arrays for colors ( shades of gray )  */
  
  unsigned char *data;
  data = malloc( iLength * sizeof(unsigned char) );
  
  if (data == NULL )
    {
      fprintf(stderr," Could not allocate memory");
      return 1;
    }
  

  
  
  
   
  Cx = -1.36809742955000002314; 
  Cy = 0.0; 

   
  NrOfCrOrbitPoints = pow(10,);
 
  for (p=2;p<=2*period+1 ;p++)
   { 
  //printf(" clear the arrays \n");       
  ClearArray( data );
  /* ------------------ draw ----------------------- */
  //printf(" Draw %llu points of critical orbit to array \n",NrOfCrOrbitPoints);       
  diffs[p] =  DrawLimitSet(NrOfCrOrbitPoints, data, p);
  
 
 
  /* ---------- file  -------------------------------------*/
  //printf(" save  data array to the pgm file \n");
  SaveArray2pgm(data, p);
  }
 
  /* --------------free memory ---------------------*/
  free(data);

 for (p=2;p<=2*period+1 ;p++) printf(" p = %u ; aver diff = %f \n", p, diffs[p]);
 
 
 
 
 
 
  return 0;
}


