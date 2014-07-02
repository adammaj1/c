
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
  2. first it creates Fatou array which is used to store rgb color values of pixels,
  fills tha array with Fatou and after that writes the Fatou from array to pgm file.
  It alows free ( non sequential) acces to "pixels"
    
  -------------------------------------------
  Adam Majewski   fraktal.republika.pl 
 
  Sobel filter 
  Gh = sum of six values ( 3 values of matrix are equal to 0 ). Each value is = pixel_color * filter_coefficients 
    
  ----------------------------------------------
  gcc i.c -lm -Wall -o2
  gcc i.c -lm -Wall -march=native -fopenmp
  ./a.out


 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <omp.h>


const int iMax =20000;
int iMaxN;

// integer ( screen) coordinate of virtual 2D array 
// Indexes of array starts from 0 not 1 
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
// c parameter of the function fc(z)
double Cx;
double Cy; 
double complex c;
// parabolic fixed point
double complex ZA;  /* atractor ZA = ZAx + ZAy*i */
double ZAx,ZAy;
// critical point
double complex ZCr = 0.0; 
double ZCrx = 0.0;
double ZCry = 0.0; 



double AR ; //=0.05  /* PixelWidth*1.5   radius of circle around attractor ZA = target set for attracting points */
double  AR2 ; //AR*AR
double ARn; // normalized see : GiveNormalizedParameters ; normalization : phi(critical point ) = 0
double ARn2;

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
  //unsigned char iComponentL = 100;
  //unsigned char iComponentR = 80;
  unsigned char iBasin = 120;

  unsigned char iLostPixels = 100;
  unsigned char iPetalL=150;
  unsigned char iPetalR= 140;
  unsigned char iInterior = 180;
  unsigned char iJulia = 0;

 /* dynamic 1D arrays for colors ( shades of gray ) */
  unsigned char *Fatou, *Julia, *Basin, *boundary;

// ====================  functions ========================================================================


// complex quadaratic polynomial 
double complex f(double complex z, double complex c)
{
  return ( z*z +c);

}


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
unsigned int g(unsigned int _iX, unsigned int _iY)
{return (_iX + (iYmax-_iY-1)*iXmax );}


// compute color of the pixel z 
// uses global variables : t, ZA, c, ER2, AR2, iPeriodChild
unsigned char GiveFatouColor(double Zx0, double Zy0)
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


 
/*
 normalization : 
   phi(critical point ) = 0
    Attraction time(zcr) = k
   attraction time of (zcr+pixelwidth) = k+1
  so :
fix iMax
 compute k
 compute distance between (zcr+pixelwidth) and alpha after k iterations 
 change AR : d( zcr, alfa, k) < newAR < d( zcr+pixelWidth,alfa, k)

then zcr will be on boundary of level sets 

 here z_cr = 0.0
 changes iterationMax and AR 
*/
double GiveNormalizedParameters()
{
 
 int i,p;
 
// initiall point of iteration is a critical point of f
 double Zx=0.0;
 double  Zy=0.0;
 double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
 /* distance from z to Alpha  after iMax iterations */
 double dcr;
 double d1;
 
 // compute dcr
 Zx = ZCrx; // critical point z=0
 Zy = ZCry;
 Zx2=Zx*Zx;
 Zy2=Zy*Zy;
 for ( i=0; i<iMax; ++i)
  {   //  forward iteration of z under fc(z) = z^2 + c
      for(p=0; p< iPeriodChild ;++p) // iMax = period !!!!
       { 
        // not needed here if (Zx2+Zy2 > ER2) return iExterior; // bailout test 
        Zy=2*Zx*Zy + Cy;
        Zx=Zx2-Zy2 +Cx;
        Zx2=Zx*Zx;
        Zy2=Zy*Zy;
        
       }
  }
  dcr = fabs(Zx- ZAx);

 // compute d1 // next pixel from critical point 
 Zx = ZCrx; 
 Zy = ZCry+PixelWidth; // move up or dpwn , because next right/left  point has the same time  
 Zx2=Zx*Zx;
 Zy2=Zy*Zy;
 for ( i=0; i<iMax; ++i)
  {   //  forward iteration of z under fc(z) = z^2 + c
      for(p=0; p< iPeriodChild ;++p) // iMax = period !!!!
       { 
        // not needed here if (Zx2+Zy2 > ER2) return iExterior; // bailout test 
        Zy=2*Zx*Zy + Cy;
        Zx=Zx2-Zy2 +Cx;
        Zx2=Zx*Zx;
        Zy2=Zy*Zy;
         
      }
      
      
  }
  d1 = fabs(Zx- ZAx);
  printf(" dcr = %.20f < d1 = %.20f < AR = %.20f \n", dcr, d1, AR);

  // now : dcr<d1<AR
  // find ARn = normalised AR  
  ARn = dcr + (d1-dcr)/2.0;
 // now dcr<AR<d1
 printf(" dcr = %.20f < ARn = %.20f < d1 = %.20f \n", dcr, ARn, d1);
  
 return ARn; // 
 
}
 

unsigned char GivePetalLevelsColor(unsigned int iX, unsigned int iY)
{
  int i; // iteration
  
  int p; // one turn around fixed point 
  double Zx, Zy; // z = Zx+Zy*I
  double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
  
  double d2, dX, dY; /* distance from z to Alpha  */

 
  
  // from screen ( intteger ) to world ( double) coordinate
  Zy=ZyMin + iY*PixelHeight; /* map from screen to world coordinate  */
  if (fabs(Zy)<PixelHeight/2) Zy=0.0; /*  */
  Zx=ZxMin + iX*PixelWidth;
 
 
  Zx2=Zx*Zx;
  Zy2=Zy*Zy;
  dX=Zx-ZAx;
  dY=Zy-ZAy;
  d2=dX*dX+dY*dY;

  i=0; 
  while (i<iMaxN ) 
    {
      // Fall-in test : iterate z and check when point z is inside target set.
      // Note that test is done after every f^{period} not after every f !!!
      if (d2 < ARn2)  return i; //100*(i % 255); // break;
        
 
      for(p=0; p< iPeriodChild ;++p) // iMax = period !!!!
       { 
        // not needed here if (Zx2+Zy2 > ER2) return iExterior; // bailout test 
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

 return 100*(i % 2);
}



int MakePetalLevels(unsigned char GoodColor, unsigned char A[] )
{
  unsigned char color;
 unsigned int iX,iY, /* indices of 2D virtual array (image) = integer coordinate of the pixel */
    i; /* index of 1D array  */
  printf(" fill the A array \n");

  #pragma omp parallel for schedule(dynamic) private(i, iX, iY, color) shared(iXmax, iYmax)
  for(iY=0; iY<iYmax; ++iY){ 
   
    printf(" petal : row %u from %u \n",iY, iYmax);    // info 
    
    for(iX=0;iX<iXmax;++iX){ 
        i= g(iX,iY); /* compute index of 1D array from indices of 2D array , range(i) = [0,iLength-1] */
      // if (i>iLength-1) printf("i = %d error on iX = %d ; iY = %d \n",i, iX , iY); // range checking 
      if (A[i]==GoodColor) {
        color = GivePetalLevelsColor(iX, iY);
        A[i]= color; /*  */
        } 
      /*  if (Zx>0 && Zy>0) A[i]=255-A[i];    check the orientation of Z-plane by marking first quadrant */
    }
  }
return 0; 
}





int ComputeFatou(unsigned char A[])
{
 unsigned char color;
 unsigned int iX,iY, /* indices of 2D virtual array (image) = integer coordinate of the pixel */
    i; /* index of 1D array  */
 double Zx, Zy;    /* Z=Zx+Zy*i   */


//  memory image ( virtual, off-screen) 
  printf(" fill the A array \n");
/*
The use of ‘shared(variable, variable2) specifies that these variables should be shared among all the threads.
The use of ‘private(variable, variable2)’ specifies that these variables should have a seperate instance in each thread.
*/
 
#pragma omp parallel for schedule(dynamic) private(i, iX, iY, Zx, Zy, color) shared(iXmax, iYmax)

 
  for(iY=0; iY<iYmax; ++iY){ 
    Zy=ZyMin + iY*PixelHeight; /* map from screen to world coordinate  */
    if (fabs(Zy)<PixelHeight/2) Zy=0.0; /*  */
    printf(" row %u from %u \n",iY, iYmax);    // info 
    for(iX=0;iX<iXmax;++iX){ 
      Zx=ZxMin + iX*PixelWidth;
      i= g(iX,iY); /* compute index of 1D array from indices of 2D array , range(i) = [0,iLength-1] */
      // if (i>iLength-1) printf("i = %d error on iX = %d ; iY = %d \n",i, iX , iY); // range checking 
      color = GiveFatouColor(Zx, Zy);
      A[i]= color; /*  */
      /*  if (Zx>0 && Zy>0) A[i]=255-A[i];    check the orientation of Z-plane by marking first quadrant */
    }
  }

 return 0;
}

// compute boudaries in A array and save them to B  array 
int ComputeBoundariesIn(unsigned char A[], unsigned char B[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
  /* sobel filter */
  unsigned char G, Gh, Gv; 
  // boundaries are saved in B array 
 
 
 
 
  printf(" find boundaries in A array using  Sobel filter\n");   
  // #pragma omp parallel for schedule(dynamic) private(i,iY,iX,Gv,Gh,G) shared(iyMax,ixMax, ER2)
  for(iY=1;iY<iYmax-1;++iY){ 
    for(iX=1;iX<iXmax-1;++iX){ 
      Gv= A[g(iX-1,iY+1)] + 2*A[g(iX,iY+1)] + A[g(iX-1,iY+1)] - A[g(iX-1,iY-1)] - 2*A[g(iX-1,iY)] - A[g(iX+1,iY-1)];
      Gh= A[g(iX+1,iY+1)] + 2*A[g(iX+1,iY)] + A[g(iX-1,iY-1)] - A[g(iX+1,iY-1)] - 2*A[g(iX-1,iY)] - A[g(iX-1,iY-1)];
      G = sqrt(Gh*Gh + Gv*Gv);
      i= g(iX,iY); /* compute index of 1D array from indices of 2D array */
      if (G==0) {B[i]=255;} /* background */
      else {B[i]=0;}  /* boundary */
    }
  }
 
 
 
  return 0;
}


// copy boundary from B to A array 
int CopyBoundariesFromTo(unsigned char B[], unsigned char A[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
 
 
  printf("copy boundaries from Julia array to Fatou array \n");
  for(iY=1;iY<iYmax-1;++iY)
    for(iX=1;iX<iXmax-1;++iX)
      {i= g(iX,iY); if (B[i]==0) A[i]=0;}
 
 
 
  return 0;
}
 




// save "A" array to pgm file and info txt file
int SaveArray2PGMFile( unsigned char A[], double k, char* text)
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
 
 fprintf(fp,"This png image shows rectangle part of dynamic plane of discrete complex dynamical system  z(n+1) = fc(zn) \n");
 fprintf(fp," where  fc(z)= z^2 + c \n");
 fprintf(fp,"with numerical approximation of parabolic Julia set \n\n");
 
 fprintf(fp,"parameter c is a root point between iPeriodParent = %d and iPeriodOfChild  = %d hyperbolic components of Mandelbrot set \n", iPeriodParent , iPeriodChild);
 fprintf(fp,"on the end of the internal ray of parent component of Mandelbrot set with angle = 1/%d in turns \n", iPeriodChild);
 fprintf(fp," c = ( %f ; %f ) \n", Cx, Cy); 
 fprintf(fp," Image shows : %s\n\n", text);
 fprintf(fp,"parabolic alfa fixed point z = ( %f ; %f )  \n", creal(ZA), cimag(ZA)); 
 fprintf(fp," radius around parabolic fixed point AR =  %f ; Pixel width = %f   \n", AR, PixelWidth);
 fprintf(fp," iMaxN =  %d ; ARn = %f \n",iMaxN, ARn);    // info 
 fprintf(fp," EscapeRadius ER =  %f ; \n", EscapeRadius);
 fprintf(fp," Maxima number of iterations : iMax = %d \n\n", iMax);
 fprintf(fp," dynamic plane : \n");
 fprintf(fp," from ZxMin = %f to ZxMax = %f\n", ZxMin, ZxMax);
 fprintf(fp," from ZyMin = %f to ZyMax = %f\n\n", ZyMin, ZyMax);
 fprintf(fp," iWidth = %d and iHeight = %d\n", iWidth, iHeight);
 fprintf(fp," distortion = %f ; It should be 1.0 \n", ((ZxMax-ZxMin)/(ZyMax -ZyMin)/(iWidth/iHeight)));
 if (t==0) fprintf(fp," No lost pixels : t = %d ; Parameters iMax and AR are good \n",t);
    else fprintf(fp," Error : there are lost pixels : t > 0 ; you should increase iMax or make bigger AR \n\n"); 
 fprintf(fp," computations made with double type numbers \n\n");
 fprintf(fp,"use (k+AR) for file names where k is a number of file ; range(k) =[0,k] so there are (k+1) png files nad (k+1) text files \n\n");
 fprintf(fp,"made with c console program  \n");
 fprintf(fp,"Adam Majewski  \n");

 printf("File %s saved. \n", filename);
 fclose(fp); 


  return 0;
}
 


int FillContour(complex double seed,  unsigned char color, unsigned char A[])
{ 
  /* 
     fills contour with black border ( color = iJulia)  using seed point inside contour 
     and horizontal lines 
     it starts from seed point, saves max right( iXmaxLocal) and max left ( iXminLocal) interior points of horizontal line,
     in new line ( iY+1 or iY-1) it computes new interior point  : iXmidLocal=iXminLocal + (iXmaxLocal-iXminLocal)/2;
     result is stored in A array : 1D array of 1-bit colors ( shades of gray)
     it does not check if index of A array is good  so memory error is possible 
  */
 
  int iXseed = (int)((creal(seed) - ZxMin)/PixelWidth);
  int iYseed = (int)((cimag(seed) - ZyMin)/PixelHeight);
 
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
    i =g(iX,iY); /* index of A array */;
 
    /* move to right */
    while (A[i]>0 && A[i]!=iExterior) 
      { A[i]=color;
	iX+=1; 
	i=g(iX,iY);  
      }
    iXmaxLocal=iX-1;
 
    /* move to left */
    iX=iXmidLocal-1; 
    i=g(iX,iY);
    while (A[i]>0 && A[i]!=iExterior) 
      { A[i]=color;
	iX-=1; 
	i=g(iX,iY); 
      }
    iXminLocal=iX+1; 
 
 
    iY+=1; /* move up */
    iXmidLocal=iXminLocal + (iXmaxLocal-iXminLocal)/2; /* new iX inside contour */
    i=g(iXmidLocal,iY); /* index of A array */;
    if ( A[i]==iJulia)  break; /*  it should not cross the border */
 
  } while  (iY<iYmax); 
 
 
  /* ------  move down ----------------- */
  iXmidLocal=iXseed;
  iY=iYseed-1;
 
 
  do{
    iX=iXmidLocal;
    i =g(iX,iY); /* index of A array */;
 
    /* move to right */
    while (A[i]>0 && A[i]!=iExterior)  /*  */
      { A[i]=color;
	iX+=1;
	i=g(iX,iY);  
      }
    iXmaxLocal=iX-1;
 
    /* move to left */
    iX=iXmidLocal-1; 
    i=g(iX,iY);
    while (A[i]>0 && A[i]!=iExterior)  /*  */
      { A[i]=color;
	iX-=1; /* move to right */
	i=g(iX,iY);  
      }
    iXminLocal=iX+1; 
 
    iY-=1; /* move down */
    iXmidLocal=iXminLocal + (iXmaxLocal-iXminLocal)/2; /* new iX inside contour */
    i=g(iXmidLocal,iY); /* index of A array */;
    if ( A[i]==iJulia)  break; /*  it should not cross the border */
  } while  (0<iY); 
// end of horizontal fill 

 // vertical fill  
 iX = iXseed;
 iY = iYseed;
 i =g(iX,iY); /* index of A array */; 

while ( A[i]>0) {

 
 	do // up 
 	{
   		A[i]=color;
   		iY+=1; 
   		i =g(iX,iY); /* index of A array */; 
 	} while ( A[i]>0);

  iX+=1;
  iY = iYseed;	 
  i =g(iX,iY); /* index of A array */; 

  } 

 // vertical fill up 
 iX = iXseed;
 iY = iYseed;
 i =g(iX,iY); /* index of A array */; 

while ( A[i]>0) {

 
 	do // up 
 	{
   		A[i]=color;
   		iY+=1; // up 
   		i =g(iX,iY); /* index of A array */; 
 	} while ( A[i]>0);

  iX+=1; // up and to the right 
  iY = iYseed;	 
  i =g(iX,iY); /* index of A array */ 

  } 

  iX = iXseed;
 iY = iYseed;
 i =g(iX,iY); /* index of A array */; 

while ( A[i]>0) {

 
 	do // up 
 	{
   		A[i]=color;
   		iY+=1; //up 
   		i =g(iX,iY); /* index of A array */; 
 	} while ( A[i]>0);

  iX-=1; // up and to the left
  iY = iYseed;	 
  i =g(iX,iY); /* index of A array */ 

  } 


 // vertical fill down
 iX = iXseed;
 iY = iYseed;
 i =g(iX,iY); /* index of A array */; 

while ( A[i]>0) {

 
 	do //  
 	{
   		A[i]=color;
   		iY-=1; // down 
   		i =g(iX,iY); /* index of A array */; 
 	} while ( A[i]>0);

  iX+=1; // down and to the right 
  iY = iYseed;	 
  i =g(iX,iY); /* index of A array */ 

  } 

  iX = iXseed;
 iY = iYseed;
 i =g(iX,iY); /* index of A array */; 

while ( A[i]>0) {

 
 	do //  
 	{
   		A[i]=color;
   		iY-=1; // down 
   		i =g(iX,iY); /* index of A array */; 
 	} while ( A[i]>0);

  iX-=1; // down and to the left
  iY = iYseed;	 
  i =g(iX,iY); /* index of A array */ 

  } 



  return 0;
}

// immediate basin of parabolic fixed point :
// fixed point is a root point of period components
//  flower has period components 
// centers of components are ZCr and it's period images 
// algorithm : fill period contours in array A with color 
int MarkImmediateBasin( int period, unsigned char color,  unsigned char A[])
{
  double complex  Z;
  int p; 
  

  // first center is the critical point
  Z = ZCr;
  
    
  for (p=0; p<period; ++p)
  {
    Z = f(Z,c); // next center is the image of critical point under f function 
    FillContour(Z , color, A); // fill next 
  }
  
  return 0; 

}
 

int setup()
{

iLength = iWidth*iHeight;/* length of array in bytes = number of bytes = number of pixels of image * number of bytes of color */
 /* dynamic 1D arrays for colors ( shades of gray ) */
  Fatou = malloc( iLength * sizeof(unsigned char) );
  Julia = malloc( iLength * sizeof(unsigned char) );
  boundary = malloc( iLength * sizeof(unsigned char) );
  Basin = malloc( iLength * sizeof(unsigned char) );
  if (Fatou == NULL || Julia==NULL || boundary==NULL || Basin==NULL)
    {
      fprintf(stderr," Could not allocate memory");
      getchar(); 
      return 1;
    }
  else printf(" memory is OK\n");


  PixelWidth=(ZxMax-ZxMin)/ (iWidth-1);
  PixelHeight=(ZyMax-ZyMin)/(iHeight-1);
  
  
 

  denominator = iPeriodChild;
  InternalAngle = 1.0/((double) denominator);
  // find the c parameter with specyfic features
  c = GiveC(InternalAngle, 1.0, iPeriodParent) ; // internal radius= 1.0 gives root point = parabolic parameter   
  Cx=creal(c);
  Cy=cimag(c);
 
  
 ZA = GiveAlfaFixedPoint( c);
 ZAx=creal(ZA);
 ZAy = cimag(ZA);

 ER2=EscapeRadius*EscapeRadius;
 AR = PixelWidth*5.0; /*   radius of circle around attractor ZA = target set for attracting points */
 AR2  = AR*AR;
 ARn = GiveNormalizedParameters(); // compute normalized AR = ARn 
 ARn2 =ARn*ARn;
 iMaxN = 10000*iMax;

 iXmax = iWidth; // range(iX) = [0, iXmax-1] 
 iYmax = iHeight; // range(iY) = [0, iYmax -1]




return 0; 
}

/* --------------------------------------------------------------------------------------------------------- */

int main(){
  
  
 printf(" setup \n");
 setup();

 // use (k+AR) for file names where k is a number of file ; range(k) =[0,k] so there are (k+1) png files nad (k+1) text files 
 printf(" mark Fatou components \n");
 ComputeFatou(Fatou);
 SaveArray2PGMFile(Fatou, 0.0+AR," Fatou components  ");

 printf(" find the Julia set as the boundary betweeen Fatou components\n");
 ComputeBoundariesIn(Fatou, Julia);
 SaveArray2PGMFile(Julia, 1.0+AR," Julia set ");

 

 printf(" Mark the immediate basin of parabolic fixed point alfa \n");
 MarkImmediateBasin(iPeriodChild, iBasin, Fatou );
 SaveArray2PGMFile(Fatou, 4.0+AR, "Fatou components, Julia set, immediate basin of parabolic fixed point alfa ");

 printf(" attraction time to parabolic fixed point alfa in an attracting petal\n");
 MakePetalLevels(iBasin , Fatou );
 SaveArray2PGMFile(Fatou, 5.0+AR, "Fatou components, Julia set, immediate basin of parabolic fixed point alfa, attraction time to parabolic fixed point alfa in one attracting petal ");

 //printf(" copy the Julia set from the Julia array to the Fatou array\n");
 //CopyBoundariesFromTo(Julia, Fatou);         
 //CopyBoundariesFromTo(Julia, Basin);
 //SaveArray2PGMFile(Fatou,  2.0+AR," Fatou components and Julia set ");
 //SaveArray2PGMFile(Basin, 3.0+AR," immediate Basin");

 printf(" find the Julia set as the boundary betweeen Fatou components\n");
 ComputeBoundariesIn(Fatou, boundary);
 SaveArray2PGMFile(boundary, 6.0+AR," boundaries ");


 //ComputeBoundariesIn(Basin, boundary);
 //printf(" copy boundaries of level sets from the Julia array to the Fatou array\n");
 //CopyBoundariesFromTo(boundary, Fatou);

 //printf("  "); 
 //SaveArray2PGMFile(Fatou, 6.0+AR," Fatou components, Julia set, immediate basin of parabolic fixed point alfa, attraction time to parabolic fixed point alfa in one attracting petal with level sets ");
 //SaveArray2PGMFile(Julia, 7.0+AR," Fatou components, Julia set, immediate basin of parabolic fixed point alfa, attraction time to parabolic fixed point alfa in one attracting petal with level sets ");
 //SaveArray2PGMFile(boundary, 8.0+AR," boundary");
 
 printf("free memory \n");
 free(Fatou);
 free(Julia);
 free(boundary); 

 return 0;

}
