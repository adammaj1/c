/*

Adam Majewski
fraktal.republika.pl

c console progam using 
* symmetry
* openMP

It uses modified DEM method (different from Milnor's) to draw parabolic julia sets



gcc t.c -lm -Wall -fopenmp -march=native 
time ./a.out


File s1000000f1.5.pgm saved. 
Cx  = 0.356763 
Cy  = 0.328582 
alfax  = 0.154508 
alfay  = 0.475528 
distorsion of image  = 1.000000 

real	123m19.106s
----------------------------
File s10000000f1.5.pgm saved. 
Cx  = 0.356763 
Cy  = 0.328582 
alfax  = 0.154508 
alfay  = 0.475528 
distorsion of image  = 1.000000 

real	1023m48.596s = 17 fours
-------------------------------------
so s100000000f1.5.pgm shoud take 17 days !!!!!! 






*/



#include <stdio.h>
#include <stdlib.h> // malloc
#include <string.h> // strcat
#include <math.h> // M_PI; needs -lm also 
#include <complex.h>
#include <omp.h> // OpenMP; needs also -fopenmp


/* --------------------------------- global variables and consts ------------------------------------------------------------ */


// virtual 2D array and integer ( screen) coordinate
// Indexes of array starts from 0 not 1 
unsigned int ix, iy; // var
unsigned int ixMin = 0; // Indexes of array starts from 0 not 1
unsigned int ixMax ; //
unsigned int iWidth ; // horizontal dimension of array
unsigned int ixAxisOfSymmetry  ; // 
unsigned int iyMin = 0; // Indexes of array starts from 0 not 1
unsigned int iyMax ; //
unsigned int iyAxisOfSymmetry  ; // 
unsigned int iyAbove ; // var, measured from 1 to (iyAboveAxisLength -1)
unsigned int iyAboveMin = 1 ; //
unsigned int iyAboveMax ; //
unsigned int iyAboveAxisLength ; //
unsigned int iyBelowAxisLength ; //
unsigned int iHeight = 1000; //  odd number !!!!!! = (iyMax -iyMin + 1) = iyAboveAxisLength + iyBelowAxisLength +1
// The size of array has to be a positive constant integer 
unsigned int iSize ; // = iWidth*iHeight; 


// memmory 1D array 
unsigned char *data;
// unsigned int i; // var = index of 1D array
unsigned int iMin = 0; // Indexes of array starts from 0 not 1
unsigned int iMax ; // = i2Dsize-1  = 
// The size of array has to be a positive constant integer 
// unsigned int i1Dsize ; // = i2Dsize  = (iMax -iMin + 1) =  ;  1D array with the same size as 2D array


/* world ( double) coordinate = dynamic plane */
  const double ZxMin=-1.5;
  const double ZxMax=1.5;
  const double ZyMin=-1.5;
  const double ZyMax=1.5;
  double PixelWidth; // =(ZxMax-ZxMin)/iXmax;
  double PixelHeight; // =(ZyMax-ZyMin)/iYmax;
  double distanceMax;
  double ratio ;
  double lambda=1.5;



// complex numbers of parametr plane 
double Cx; // c =Cx +Cy * i
double Cy;
double complex c; // 

double complex alfa; // alfa fixed point alfa=f(alfa)

unsigned long int iterMax  = 10000000; //iHeight*100;
double ER = 2.0; // Escape Radius for bailout test 
double ER2;


/* ------------------------------------------ functions -------------------------------------------------------------*/



/* find c in component of Mandelbrot set 
 
 uses code by Wolf Jung from program Mandel
 see function mndlbrot::bifurcate from mandelbrot.cpp
 http://www.mndynamics.com/indexp.html

  */
double complex GiveC(double InternalAngleInTurns, double InternalRadius, unsigned int period)
{
  //0 <= InternalRay<= 1
  //0 <= InternalAngleInTurns <=1
  double t = InternalAngleInTurns *2*M_PI; // from turns to radians
  double R2 = InternalRadius * InternalRadius;
  //double Cx, Cy; /* C = Cx+Cy*i */
  switch ( period ) // of component 
  {
    case 1: // main cardioid
      Cx = (cos(t)*InternalRadius)/2-(cos(2*t)*R2)/4; 
      Cy = (sin(t)*InternalRadius)/2-(sin(2*t)*R2)/4; 
      break;
   case 2: // only one component 
      Cx = InternalRadius * 0.25*cos(t) - 1.0;
      Cy = InternalRadius * 0.25*sin(t); 
      break;
  // for each period  there are 2^(period-1) roots. 
  default: // higher periods : to do
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


double GiveDistanceBetween(double complex z1, double complex z2 )
{double dx,dy;
 
 dx = creal(z1) - creal(z2);
 dy = cimag(z1) - cimag(z2);
 return sqrt(dx*dx+dy*dy);
 
} 



int setup()
{

  unsigned int period = 5; // period of secondary component joined by root point
  unsigned int denominator;
  double InternalAngle;
  

  denominator = period;
  InternalAngle = 1.0/((double) denominator);

  c = GiveC(InternalAngle, 1.0, 1) ;
  Cx=creal(c);
  Cy=cimag(c);
  alfa = GiveAlfaFixedPoint(c);

  /* 2D array ranges */
  if (!(iHeight % 2)) iHeight+=1; // it sholud be even number (variable % 2) or (variable & 1)
  iWidth = iHeight;
  iSize = iWidth*iHeight; // size = number of points in array 
  // iy
  iyMax = iHeight - 1 ; // Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].
  iyAboveAxisLength = (iHeight -1)/2;
  iyAboveMax = iyAboveAxisLength ; 
  iyBelowAxisLength = iyAboveAxisLength; // the same 
  iyAxisOfSymmetry = iyMin + iyBelowAxisLength ; 
  // ix
  
  ixMax = iWidth - 1;

  /* 1D array ranges */
  // i1Dsize = i2Dsize; // 1D array with the same size as 2D array
  iMax = iSize-1; // Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].


 /* Pixel sizes */
  PixelWidth = (ZxMax-ZxMin)/ixMax; //  ixMax = (iWidth-1)  step between pixels in world coordinate 
  PixelHeight = (ZyMax-ZyMin)/iyMax;
  ratio = ((ZxMax-ZxMin)/(ZyMax-ZyMin))/((float)iWidth/(float)iHeight); // it should be 1.000 ...
  distanceMax = PixelWidth; 
  //
  
  ER2 = ER * ER;

  /* create dynamic 1D arrays for colors ( shades of gray ) */
  
  data = malloc( iSize * sizeof(unsigned char) );
  if (data == NULL )
    {
      fprintf(stderr," Could not allocate memory\n");
      getchar();
      return 1;
    }
  else fprintf(stderr," memory is OK \n");

  
 
  return 0;

}



// from screen to world coordinate ; linear mapping
// uses global cons
double GiveZx(unsigned int ix)
{ return (ZxMin + ix*PixelWidth );}

// uses globaal cons
double GiveZy(unsigned int iy)
{ return (ZyMax - iy*PixelHeight);} // reverse y axis


// forward iteration of complex quadratic polynomial
// fc(z) = z*z +c
// z0 = critical point 
// uses global var 
long int GiveLastIteration(double Zx, double Zy )
{ 
    long int iter;
    // z= Zx + Zy*i
    double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
 
     /* initial value of orbit = critical point Z= 0 */
                       
                        Zx2=Zx*Zx;
                        Zy2=Zy*Zy;
                        // for iter from 0 to iterMax because of ++ after last loop
                        for (iter=0; iter<iterMax && ((Zx2+Zy2)<ER2); ++iter)
                        {
                            Zy=2*Zx*Zy + Cy;
                            Zx=Zx2-Zy2 + Cx;
                            Zx2=Zx*Zx;
                            Zy2=Zy*Zy;
                        };
  
   return iter;
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
-----------------


"
The algorithm I am using is a mix of techinques I have been taught by other people.

This picture was drawn by scanning all the pixels and computing the color of each of them by iterating the corresponding complex number and seeing what happens.

One of the technique is very general in its applications. Let P be the function which you want to computer the Julia set, when you inductively compute
z_n+1=P(z_n)=P^n(z_0),
also compute the derivative
d_n=(P^n)'(z_0)
It can be done inductively by the chain rule: d_n+1=P'(z_n).d_n

Then the idea is to pretend that the image of the pixel of size epsilon by P^n is roughly of size epsilon.d_n. This is not always accurate but works surprisingly well in the pictures.

Then the trick is to choose a point 'p' that we know to belong to the Julia set. Now when you iterate, test wether the distance from z_n to 'p' is less than epsilon.d_n If this is so, then there is a good chance that z_0 is within distance of order epsilon to the Julia set, so color your point in black (or whichever color you chose for the Julia set).

Usually, 'p'= a repelling fixed point of P is a good choice.

For Siegel disks and Herman rings, a good choice of 'p' is the critical point 'c'.

Namely,
if |z_n-c|<epsilon*d_n
then colored z_0 in black.
 
Your Herman ring can be enhanced a lot using the tricks which described above."

Regards,
Yang Fei

A pixel of size epsilon, after iterating n times by P, 
we obtain a domain is roughly of size epsilon*d_n, where d_n means the derivative of P^n at z_0.


-----------------------------
"The following method was taught to me by Christian Henriksen.
It is a distance estimator method (different from Milnor's).
When iterating forward, keep also tract of the value of the derivative.
More precisely when computing z_n from z_0, also compute dz_n/dz_0.
by the chain rule, it is equal to the product of the derivatives of f taken at z_0, z_1, .. , up to z_n-1
Therefore you can compute it recursively alongside and like the value of z_n.
Then choose a point in the Julia set, here I chose the critical point.
Then if the ball :
- of center z_n 
- and radius lambda * size of a pixel * derivative 
contains this point there is a good chance 
that the pixel of center z_0 
contains an iterated preimage of the critical point, thus I color z_0 in black.

Here lambda is a thickness factor that you must choose. 
not too small otherwise you do not see enough points, but not too big because otherwise you get artifacts.

Best Regards,
Arnaud."


 */
 int jdist(double Zx, double Zy, double Cx, double Cy ,  int iter_max)
 { 
 int i;
 double x = Zx, /* Z = x+y*i */
         y = Zy, 
         /* Zp = xp+yp*1 = 1  */
         xp = 1, 
         yp = 0, 
         /* temporary */
         nz,  
         nzp,
         /* a = abs(z) */
         //a,
         radius,
         distance; 


 
 for (i = 1; i <= iter_max; i++)
  { /* first derivative   zp = 2*z*zp  = xp + yp*i; */
    
    nz = 2*(x*xp - y*yp) ; 
    yp = 2*(x*yp + y*xp); 
    xp = nz;
    /* z = z*z + c = x+y*i */
    nz = x*x - y*y + Cx; 
    y = 2*x*y + Cy; 
    x = nz; 
    //
    nz = (x*x + y*y); // abs2(z)
    nzp = (xp*xp + yp*yp);
    if ( nz > ER2) return 0; // if escapes //nzp > 1e60 ||
    if (nzp>1e60) return 0; // ? interior ?
    distance=GiveDistanceBetween(alfa,x+y*I);
    nzp = sqrt(xp*xp + yp*yp);
    radius=lambda*PixelWidth*nzp;
    
    if (distance<radius)return 1; // 
  }
 
 return 0; // 
 }


unsigned char GiveColor(unsigned int ix, unsigned int iy)
{ 
   double Zx, Zy; //  Z= Zx+ZY*i;
  unsigned char color; // gray from 0 to 255 
// unsigned char ColorList[]={255,230,180}; /* shades of gray used in imaage
  //color=0;
  // from screen to world coordinate 
  Zx = GiveZx(ix);
  Zy = GiveZy(iy);
  //if ( GiveLastIteration(Zx, Zy ) == iterMax)
    // color = 0; // interior
    // else { 

// only dist
  if ( GiveDistanceBetween(alfa,Zx+Zy*I)<0.3) lambda=3.0; else lambda=1.5;
  if (jdist(Zx, Zy, Cx, Cy , iterMax))
                 color = 0; // black = boundary = Julia set
                 else color = 255; // Fatou set
          //}
  return color;
}


/* -----------  array functions -------------- */


/* gives position of 2D point (iX,iY) in 1D array  ; uses also global variable iWidth */
unsigned int Give_i(unsigned int ix, unsigned int iy)
{ return ix + iy*iWidth; }
//  ix = i % iWidth;
//  iy = (i- ix) / iWidth;
//  i  = Give_i(ix, iy);




// plots raster point (ix,iy) 
int PlotPoint(unsigned int ix, unsigned int iy, unsigned char iColor)
{
 unsigned i; /* index of 1D array */
 i = Give_i(ix,iy); /* compute index of 1D array from indices of 2D array */
 data[i] = iColor;

return 0;
}


// fill array 
// uses global var :  ...
// scanning complex plane 
int FillArray(unsigned char data[] )
{
  unsigned int ix, iy; // pixel coordinate 


// for all pixels of image 
for(iy = iyMin; iy<=iyMax; ++iy) 
  { printf(" %d z %d\n", iy, iyMax); //info 
    for(ix= ixMin; ix<=ixMax; ++ix) PlotPoint(ix, iy, GiveColor(ix, iy) ); //  
   } 
   
 return 0;
}


// fill array using symmetry of image 
// uses global var :  ...
int FillArraySymmetric(unsigned char data[] )
{
   
 unsigned char Color; // gray from 0 to 255 

printf("axis of symmetry \n"); 
iy = iyAxisOfSymmetry; 
#pragma omp parallel for schedule(dynamic) private(ix,Color) shared(ixMin,ixMax, iyAxisOfSymmetry)
for(ix=ixMin;ix<=ixMax;++ix) {//printf(" %d from %d\n", ix, ixMax); //info  
                              PlotPoint(ix, iy, GiveColor(ix, iy));
}


/*
The use of ‘shared(variable, variable2) specifies that these variables should be shared among all the threads.
The use of ‘private(variable, variable2)’ specifies that these variables should have a seperate instance in each thread.
*/

#pragma omp parallel for schedule(dynamic) private(iyAbove,ix,iy,Color) shared(iyAboveMin, iyAboveMax,ixMin,ixMax, iyAxisOfSymmetry)

// above and below axis 
for(iyAbove = iyAboveMin; iyAbove<=iyAboveMax; ++iyAbove) 
  {printf(" %d from %d\n", iyAbove, iyAboveMax); //info 
  for(ix=ixMin; ix<=ixMax; ++ix) 

  { // above axis compute color and save it to the array
    iy = iyAxisOfSymmetry + iyAbove;
    Color = GiveColor(ix, iy);
    PlotPoint(ix, iy, Color ); 
    // below the axis only copy Color the same as above without computing it 
    PlotPoint(ixMax-ix, iyAxisOfSymmetry - iyAbove , Color ); 
   } 
}  
 return 0;
}




// Check Orientation of image : first quadrant in upper right position
// uses global var :  ...
int CheckOrientation(unsigned char data[] )
{
   unsigned int ix, iy; // pixel coordinate 
   double Zx, Zy; //  Z= Zx+ZY*i;
   unsigned i; /* index of 1D array */
   for(iy=iyMin;iy<=iyMax;++iy) 
   {
     Zy = GiveZy(iy);
    for(ix=ixMin;ix<=ixMax;++ix) 
   {

    // from screen to world coordinate 
    Zx = GiveZx(ix);
     i = Give_i(ix, iy); /* compute index of 1D array from indices of 2D array */
     if (Zx>0 && Zy>0) data[i]=255-data[i];   // check the orientation of Z-plane by marking first quadrant */

    }
   }
   
  return 0;
}



// save data array to pgm file 
int SaveArray2PGMFile( unsigned char data[])
{
  
  FILE * fp;
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  char name [10]; /* name of file */
  sprintf(name,"s%luf%3.1f", iterMax, lambda); /*  */
  char *filename =strcat(name,".pgm");
  char *comment="# ";/* comment should start with # */

  /* save image to the pgm file  */      
  fp= fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P5\n %s\n %u %u\n %u\n",comment,iWidth,iHeight,MaxColorComponentValue);  /*write header to the file*/
  fwrite(data,iSize,1,fp);  /*write image data bytes to the file in one step */
  printf("File %s saved. \n", filename);
  fclose(fp);

  return 0;
}


int info()
{
 // diplay info messages
  printf("Cx  = %f \n", Cx); 
  printf("Cy  = %f \n", Cy);
  // 
 
  printf("alfax  = %f \n", creal(alfa));
  printf("alfay  = %f \n", cimag(alfa));
  printf("distorsion of image  = %f \n", ratio);
return 0;
}


/* -----------------------------------------  main   -------------------------------------------------------------*/
int main()
{
// here are procedures for creating image file

  setup();
  // compute colors of pixels = image
  //FillArray( data ); // no symmetry
  FillArraySymmetric(data); 
  //  CheckOrientation( data );
  SaveArray2PGMFile( data); // save array (image) to pgm file 
  free(data);
  info();
  return 0;
}
