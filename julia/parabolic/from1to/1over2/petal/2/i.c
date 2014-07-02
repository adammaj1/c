
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



const int iMax =20000;

/* iXmax/iYmax =  */
//#define iXmax 2000 /* height of image in pixels */
//#define iYmax 1000


// virtual 2D array and integer ( screen) coordinate
// Indexes of array starts from 0 not 1 
//unsigned int ix, iy; // var
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



/* fc(z) = z*z + c */

#define iPeriodChild 2 // Period of secondary component joined by root point with the parent component 
int iPeriodParent = 1; // main cardioid of Mandelbrot set
unsigned int denominator;/* denominator of internal angle */
double InternalAngle;
//
double Cx;
double Cy; 
double complex c;
double complex ZA;  /* atractor ZA = ZAx + ZAy*i */
double ZAx,ZAy;



double AR ; //=0.05  /* PixelWidth*1.5   radius of circle around attractor ZA = target set for attracting points */
double  AR2 ; //AR*AR
double EscapeRadius=2.0; /* radius of circle around origin; its complement is a target set for escaping points */
double ER2; //=EscapeRadius*EscapeRadius;
//#define alfa (1-sqrt(1-4*Cx))/2 /* attracting or parabolic fixed point z = alfa */
//#define beta (1+sqrt(1-4*Cx))/2 /* repelling or parabolic fixed point z = beta */
unsigned int t=0; // number of Lost Pixels ; it should be zero ! 


/* color */
  
 /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  unsigned char iExterior = 245;
  unsigned char iInteriorL = 200;
  unsigned char iInteriorR = 220;
  unsigned char iComponentL = 100;
  unsigned char iComponentR = 80;

  unsigned char iLostPixels = 100;
  unsigned char iPetalL=150;
  unsigned char iPetalR= 140;
  unsigned char iInterior = 180;
  unsigned char iJulia = 0;

 /* dynamic 1D arrays for colors ( shades of gray ) */
  unsigned char *data, *edge;

// ====================  functions ========================================================================

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


// compute color of the pixel z 
// uses global variables : t, ZA, c, ER2, AR2, iPeriodChild
unsigned char GiveColor(double Zx0, double Zy0)
{
  int i; // iteration
  int p; // one turn around fixed point 
  double Zx, Zy; // z = Zx+Zy*I
  double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
  double d2, dX, dY; /* distance from z to Alpha  */


  Zx=Zx0; /* initial value of orbit  */
  Zy=Zy0;
  Zx2=Zx*Zx;
  Zy2=Zy*Zy;
  dX=Zx-ZAx;
  dY=Zy-ZAy;
  d2=dX*dX+dY*dY;

  i=0; 
  while (i<iMax ) 
    {
      // Fall-in test : iterate z and check when point z is inside target set.
      // Note that test is done after every f^{period} not after every f !!!
      if (d2 < AR2) 
        {if (Zx>ZAx) return iInteriorR; else return iInteriorL; 
         // sd2=sqrt(d2); 
          //a= (zy/sd2);
          //a2=a*a;
          //b:(zx/sd2)+1.0;
          //b2=b*b;   
          //d=sqrt(a2+b2);
          
        }
 
      for(p=0; p< iPeriodChild ;++p) // iMax = period !!!!
       { 
        if (Zx2+Zy2 > ER2) return iExterior; // bailout test 
        Zy=2*Zx*Zy + Cy;
        Zx=Zx2-Zy2 +Cx;
        Zx2=Zx*Zx;
        Zy2=Zy*Zy;
        i+=1; 
       }
      

      
      //
      dX=Zx-ZAx;
      dY=Zy-ZAy;
      d2=dX*dX+dY*dY;
    };

 // it shouldn't happened !!!!
 t+=1; // number of iLostPixels = not escaping and not attracted points; it shoould be zero !!!!

return iLostPixels; // iMax too low ??? 
}


int FillContour(double dXseed, double dYseed,  unsigned char color, unsigned char A[])
{ 
  /* 
     fills contour with black border ( color = iJulia)  using seed point inside contour 
     and horizontal lines 
     it starts from seed point, saves max right( iXmaxLocal) and max left ( iXminLocal) interior points of horizontal line,
     in new line ( iY+1 or iY-1) it computes new interior point  : iXmidLocal=iXminLocal + (iXmaxLocal-iXminLocal)/2;
     result is stored in A array : 1D array of 1-bit colors ( shades of gray)
     it does not check if index of A array is good  so memory error is possible 
  */
 
  int iXseed = (int)((dXseed - ZxMin)/PixelWidth);
  int iYseed = (int)((dYseed - ZyMin)/PixelHeight);
 
  int iX, /* seed integer coordinate */
    iY=iYseed,
    /* most interior point of line iY */
    iXmidLocal=iXseed, 
    /* min and max of interior points of horizontal line iY */
    iXminLocal, 
    iXmaxLocal; 
  int i ; /* index of A array */;
 
 // fill horizontal lines , first up from iXseed, later down from iXseed
  /* ---------  move up --------------- */ 
  do{
    iX=iXmidLocal;
    i =f(iX,iY); /* index of A array */;
 
    /* move to right */
    while (A[i]>0 && A[i]!=iExterior) 
      { A[i]=color;
	iX+=1; 
	i=f(iX,iY);  
      }
    iXmaxLocal=iX-1;
 
    /* move to left */
    iX=iXmidLocal-1; 
    i=f(iX,iY);
    while (A[i]>0 && A[i]!=iExterior) 
      { A[i]=color;
	iX-=1; 
	i=f(iX,iY); 
      }
    iXminLocal=iX+1; 
 
 
    iY+=1; /* move up */
    iXmidLocal=iXminLocal + (iXmaxLocal-iXminLocal)/2; /* new iX inside contour */
    i=f(iXmidLocal,iY); /* index of A array */;
    if ( A[i]==iJulia)  break; /*  it should not cross the border */
 
  } while  (iY<iYmax); 
 
 
  /* ------  move down ----------------- */
  iXmidLocal=iXseed;
  iY=iYseed-1;
 
 
  do{
    iX=iXmidLocal;
    i =f(iX,iY); /* index of A array */;
 
    /* move to right */
    while (A[i]>0 && A[i]!=iExterior)  /*  */
      { A[i]=color;
	iX+=1;
	i=f(iX,iY);  
      }
    iXmaxLocal=iX-1;
 
    /* move to left */
    iX=iXmidLocal-1; 
    i=f(iX,iY);
    while (A[i]>0 && A[i]!=iExterior)  /*  */
      { A[i]=color;
	iX-=1; /* move to right */
	i=f(iX,iY);  
      }
    iXminLocal=iX+1; 
 
    iY-=1; /* move down */
    iXmidLocal=iXminLocal + (iXmaxLocal-iXminLocal)/2; /* new iX inside contour */
    i=f(iXmidLocal,iY); /* index of A array */;
    if ( A[i]==iJulia)  break; /*  it should not cross the border */
  } while  (0<iY); 
// end of horizontal fill 

 // vertical fill  
 iX = iXseed;
 iY = iYseed;
 i =f(iX,iY); /* index of A array */; 

while ( A[i]>0) {

 
 	do // up 
 	{
   		A[i]=color;
   		iY+=1; 
   		i =f(iX,iY); /* index of A array */; 
 	} while ( A[i]>0);

  iX+=1;
  iY = iYseed;	 
  i =f(iX,iY); /* index of A array */; 

  } 

 // vertical fill up 
 iX = iXseed;
 iY = iYseed;
 i =f(iX,iY); /* index of A array */; 

while ( A[i]>0) {

 
 	do // up 
 	{
   		A[i]=color;
   		iY+=1; // up 
   		i =f(iX,iY); /* index of A array */; 
 	} while ( A[i]>0);

  iX+=1; // up and to the right 
  iY = iYseed;	 
  i =f(iX,iY); /* index of A array */ 

  } 

  iX = iXseed;
 iY = iYseed;
 i =f(iX,iY); /* index of A array */; 

while ( A[i]>0) {

 
 	do // up 
 	{
   		A[i]=color;
   		iY+=1; //up 
   		i =f(iX,iY); /* index of A array */; 
 	} while ( A[i]>0);

  iX-=1; // up and to the left
  iY = iYseed;	 
  i =f(iX,iY); /* index of A array */ 

  } 


 // vertical fill down
 iX = iXseed;
 iY = iYseed;
 i =f(iX,iY); /* index of A array */; 

while ( A[i]>0) {

 
 	do //  
 	{
   		A[i]=color;
   		iY-=1; // down 
   		i =f(iX,iY); /* index of A array */; 
 	} while ( A[i]>0);

  iX+=1; // down and to the right 
  iY = iYseed;	 
  i =f(iX,iY); /* index of A array */ 

  } 

  iX = iXseed;
 iY = iYseed;
 i =f(iX,iY); /* index of A array */; 

while ( A[i]>0) {

 
 	do //  
 	{
   		A[i]=color;
   		iY-=1; // down 
   		i =f(iX,iY); /* index of A array */; 
 	} while ( A[i]>0);

  iX-=1; // down and to the left
  iY = iYseed;	 
  i =f(iX,iY); /* index of A array */ 

  } 



  return 0;
}

int ComputeBoundariesIn(unsigned char A[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
  /* sobel filter */
  unsigned char G, Gh, Gv; 
  // boundaries are in edge array ( global var )
 
 
 
 
  printf(" find boundaries in A array using  Sobel filter\n");   
  // #pragma omp parallel for schedule(dynamic) private(i,iY,iX,Gv,Gh,G) shared(iyMax,ixMax, ER2)
  for(iY=1;iY<iyMax-1;++iY){ 
    for(iX=1;iX<ixMax-1;++iX){ 
      Gv= A[Give_i(iX-1,iY+1)] + 2*A[Give_i(iX,iY+1)] + A[Give_i(iX-1,iY+1)] - A[Give_i(iX-1,iY-1)] - 2*A[Give_i(iX-1,iY)] - A[Give_i(iX+1,iY-1)];
      Gh= A[Give_i(iX+1,iY+1)] + 2*A[Give_i(iX+1,iY)] + A[Give_i(iX-1,iY-1)] - A[Give_i(iX+1,iY-1)] - 2*A[Give_i(iX-1,iY)] - A[Give_i(iX-1,iY-1)];
      G = sqrt(Gh*Gh + Gv*Gv);
      i= Give_i(iX,iY); /* compute index of 1D array from indices of 2D array */
      if (G==0) {edge[i]=255;} /* background */
      else {edge[i]=0;}  /* boundary */
    }
  }
 
 
 
  return 0;
}
 
 
int CopyBoundariesTo(unsigned char A[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
 
 
  printf("copy boundaries from edge array to data array \n");
  for(iY=1;iY<iyMax-1;++iY)
    for(iX=1;iX<ixMax-1;++iX)
      {i= Give_i(iX,iY); if (edge[i]==0) A[i]=0;}
 
 
 
  return 0;
}
 

// save "A" array to pgm file and info txt file
int SaveArray2PGMFile( unsigned char A[], double k)
{
 
  FILE * fp;
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  char name [100]; /* name of file  , error Naruszenie ochrony pamięci  */
  sprintf(name,"%f", k); /*  */
  char *filename =strcat(name,".pgm");
  char *comment="# ";/* comment should start with # */
 
  /* save image to the pgm file  */      
  fp= fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P5\n %s\n %u %u\n %u\n",comment,iWidth,iHeight,MaxColorComponentValue);  /*write header to the file*/
  fwrite(A,iLength,1,fp);  /*write A array to the file in one step */
  printf("File %s saved. \n", filename);
  fclose(fp);

  // save info text filename
 filename =strcat(name,".txt");
 fp= fopen(filename,"wb");
 fprintf(fp,"This image shows rectangle part of dynamic plane of discrete complex dynamical system  z(n+1) = fc(zn) \n");
 fprintf(fp," where  fc(z)= z^2 + c \n");
 fprintf(fp,"with numerical approximation of parabolic Julia set \n");
 fprintf(fp,"parameter c is a root point between iPeriodParent = %d and iPeriodOfChild  = %d hyperbolic components of Mandelbrot set \n", iPeriodParent , iPeriodChild);
 fprintf(fp,"on the end of the internal ray with angle = 1/%d in turns \n", iPeriodChild);
 fprintf(fp," c = ( %f ; %f ) \n", Cx, Cy); 
 fprintf(fp,"parabolic alfa fixed point z = ( %f ; %f )  \n", creal(ZA), cimag(ZA)); 
 fprintf(fp," AR =  %f ; Pixel width = %f   \n", AR, PixelWidth);
 fprintf(fp," EscapeRadius ER =  %f ; \n", EscapeRadius);
 fprintf(fp," Maxima number of iterations : iMax = %d \n", iMax); 
 if (t==0) fprintf(fp," No lost pixels : t = %d \n",t);
    else fprintf(fp," Error : there are lost pixels : t > 0 ; you should increase iMax or make bigger AR \n"); 
 fprintf(fp," computations made with double type numbers \n");
 printf("File %s saved. \n", filename);
 fclose(fp); 


  return 0;
}
 

 

/* --------------------------------------------------------------------------------------------------------- */

int main(){
  
   unsigned char color;
   unsigned int iX,iY, /* indices of 2D virtual array (image) = integer coordinate of the pixel */
    i; /* index of 1D array  */
    
  
  /* */
  double Zx, Zy;    /* Z=Zx+Zy*i   */
  /* sobel filter */
  unsigned char G, Gh, Gv; 


  
  
 printf(" setup \n");
 iLength = iWidth*iHeight;/* length of array in bytes = number of bytes = number of pixels of image * number of bytes of color */
 /* dynamic 1D arrays for colors ( shades of gray ) */
  data = malloc( iLength * sizeof(unsigned char) );
  edge = malloc( iLength * sizeof(unsigned char) );
  if (data == NULL || edge==NULL)
    {
      fprintf(stderr," Could not allocate memory");
      getchar(); 
      return 1;
    }
  else printf(" memory is OK\n");


  
  PixelWidth=(ZxMax-ZxMin)/ (iWidth-1);
  PixelHeight=(ZyMax-ZyMin)/(iHeight-1);
  ER2=EscapeRadius*EscapeRadius;
  AR = PixelWidth*20.0; /*   radius of circle around attractor ZA = target set for attracting points */
  AR2  = AR*AR;

  denominator = iPeriodChild;
  InternalAngle = 1.0/((double) denominator);
 
  c = GiveC(InternalAngle, 1.0, iPeriodParent) ; // internal radius= 1.0 gives root point = parabolic parameter   
  Cx=creal(c);
  Cy=cimag(c);
 
  
 ZA = GiveAlfaFixedPoint( c);
 ZAx=creal(ZA);
 ZAy = cimag(ZA);

 iXmax = iWidth; // range(iX) = [0, iXmax-1] 
 iYmax = iHeight; // range(iY) = [0, iYmax -1]


 //  memory image ( virtual, off-screen) 
  printf(" fill the data array \n");
  for(iY=0; iY<iYmax; ++iY){ 
    Zy=ZyMin + iY*PixelHeight; /*  */
    if (fabs(Zy)<PixelHeight/2) Zy=0.0; /*  */
    printf(" row %u from %u \n",iY, iYmax);    // info 
    for(iX=0;iX<iXmax;++iX){ 
      Zx=ZxMin + iX*PixelWidth;
      i= f(iX,iY); /* compute index of 1D array from indices of 2D array , range(i) = [0,iLength-1] */
      // if (i>iLength-1) printf("i = %d error on iX = %d ; iY = %d \n",i, iX , iY); // range checking 
      color = GiveColor(Zx, Zy);
      data[i]= color; /*  */
      /*  if (Zx>0 && Zy>0) data[i]=255-data[i];    check the orientation of Z-plane by marking first quadrant */
    }
  }


  printf(" find boundaries in data array using  Sobel filter\n");   
  for(iY=1; iY<iYmax-1; ++iY) // not use first and last point because they have not all neighbours !!
   { 
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
    for(iY=1; iY<iYmax-1; ++iY){ 
     for(iX=1;iX<iXmax-1;++iX)
      {i= f(iX,iY); /* compute index of 1D array from indices of 2D array */
  	if (edge[i]==0) data[i]=0;}}


  printf(" save  data array to the file \n");
  SaveArray2PGMFile(data,AR);
  printf(" Mark the immediate basin of parabolic fixed point alfa \n");
  FillContour(ZAx -0.4, ZAy,  iComponentR, data);
  FillContour(ZAx+0.4, ZAy,  iComponentL, data);
  SaveArray2PGMFile(data,1.0+AR);
  /* --------------free memory ---------------------*/
  free(data);
  free(edge);
  

  return 0;
}
