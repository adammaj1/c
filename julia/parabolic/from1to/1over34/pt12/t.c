/*

  Adam Majewski
  adammaj1 aaattt o2 dot pl  // o like oxygen not 0 like zero 
  fraktal.republika.pl

  c console progam 

  gcc t.c -lm -Wall -march=native 
  time ./a.out



  method of filling gaps in rays turns 
  should be the same as for drawing rays !!!!!!!!!!!!!

  !!!!!---------- algorithm ---------!!!!!!!!
  trap ( target set) for iteration :
  * of interior points = interior of circle centered at alfa and radius = dMaxDistance2Alfa, with iPeriodChild sectors described by periodic external rays landing on alfa
  * of exterior points = exterior of circle centered at origin ( z=0) and radius = ER
  Iterate point until it fails to one of 2 traps.
  If it it fail into interior trap then check in which sector is it. Color = number of sector

  Interior target set should be all inside Julia set 
  ( of course not counting external rays that cross it and exterior parts which are very thin = means width < pixelwidth)
  if it is a circle around alfa then it shrinks as iPeriodOfChild grows
  ( = time of creating image grows)
  -----------------------------------------
  Better is to take as a interior target set :
  part of this circle bordered by 2 external rays 
  (One can choose the longest one !!!)
  and color the interior proportionaly to 
  (i % iPeriodOfChild)


 
Numerical approximation of parabolic Julia set for fc(z)= z^2 + c 
iPeriodParent = 1 
iPeriodOfChild  = 34 
parameter c = ( 0.258368 ; 0.001564 ) 
alfa fixed point z = ( 0.491487 ; 0.091875 )  
Image Width = 3.000000 
PixelWidth = 0.001500 
size of target set in screen units = iMaxDistance2Alfa  = 280 pixels 
size of target set in world units = dMaxDistance2Alfa  = 0.420000 ; 
Maximal number of iterations = iterMax = 10000000 
ratio of image  = 1.000000 ; it should be 1.000 ...

real	44m17.040s
user	44m15.250s
sys	0m0.320s

















































t = 0.0000000001
t = 0.0000000001
t = 0.0000000002
t = 0.0000000005
t = 0.0000000009
t = 0.0000000019
t = 0.0000000037
t = 0.0000000075
t = 0.0000000149
t = 0.0000000298
t = 0.0000000596
t = 0.0000001192
t = 0.0000002384
t = 0.0000004768
t = 0.0000009537
t = 0.0000019073
t = 0.0000038147
t = 0.0000076294
t = 0.0000152588
t = 0.0000305176
t = 0.0000610352
t = 0.0001220703
t = 0.0002441406
t = 0.0004882813
t = 0.0009765625
t = 0.0019531250
t = 0.0039062500
t = 0.0078125000
t = 0.0156250000
t = 0.0312500000
t = 0.0625000000
t = 0.1250000000
t = 0.2500000000
t = 0.5000000000
t = 0.0000000001

*/



#include <stdio.h>
#include <stdlib.h> // malloc
#include <string.h> // strcat
#include <math.h> // M_PI; needs -lm also 
#include <complex.h>



/* --------------------------------- global variables and consts ------------------------------------------------------------ */


#define iPeriodChild 34 // Period of secondary component joined by root point with the parent component 
#define numerator   1 // of internal angle in turns
int iPeriodParent = 1; // main cardioid of Mandelbrot set

// radius of the target set ( circle around alfa fixed point ); it is related with iHeight
// so changing iHeight needs change of iMaxDistance2Alfa
#define iMaxDistance2Alfa  200 // distance point to alfa fixed point in pixels  150 when iHeight=1000; 280 when iHeight=2000, note also iterMax and  iHeight
int iMaxDistance2Alfa2;
double dMaxDistance2Alfa2; // = (iMaxDistance2Alfa*PixelWidth)^2
double dMaxDistance2Alfa;

// virtual 2D array and integer ( screen) coordinate
// Indexes of array starts from 0 not 1 
//unsigned int ix, iy; // var
static unsigned int ixMin = 0; // Indexes of array starts from 0 not 1
static unsigned int ixMax ; //
static unsigned int iWidth ; // horizontal dimension of array

static unsigned int iyMin = 0; // Indexes of array starts from 0 not 1
static unsigned int iyMax ; //

static unsigned int iHeight = 4000; //  odd number !!!!!! = (iyMax -iyMin + 1) = iyAboveAxisLength + iyBelowAxisLength +1
// The size of array has to be a positive constant integer 
static unsigned int iSize ; // = iWidth*iHeight; 


// memmory 1D array 
unsigned char *rays;
unsigned char *data;
unsigned char *edge;

// unsigned int i; // var = index of 1D array
//static unsigned int iMin = 0; // Indexes of array starts from 0 not 1
static unsigned int iMax ; // = i2Dsize-1  = 
// The size of array has to be a positive constant integer 
// unsigned int i1Dsize ; // = i2Dsize  = (iMax -iMin + 1) =  ;  1D array with the same size as 2D array


/* world ( double) coordinate = dynamic plane */
static   const double ZxMin=-1.5;
static  const double ZxMax=1.5;
static  const double ZyMin=-1.5;
static  const double ZyMax=1.5;
static  double PixelWidth; // =(ZxMax-ZxMin)/ixMax;
static  double PixelHeight; // =(ZyMax-ZyMin)/iyMax;
static  double ratio ;
 


// complex numbers of parametr plane 
double Cx; // c =Cx +Cy * i
double Cy;
double complex c; // parameter of function fc(z)=z^2 + c

double complex alfa; // alfa fixed point alfa=f(alfa)
double dAlfaX, dAlfaY;
int iAlfaX, iAlfaY;


static unsigned long int iterMax  = 1000000; //iHeight*100;

static double ER = 2.0; // Escape Radius for bailout test 
static double ER2;


/* colors = shades of gray from 0 to 255 */
// 8 bit color = int number from 0 to 255
unsigned char iColorsOfInterior[iPeriodChild]; //={255,230,180, 160,140,120,100}; // NumberOfPetal of colors = iPeriodChild
static unsigned char iColorOfExterior = 245;
unsigned char iInterior = 200;


// array with angles in turns of points of periodic rays landing on alfa fixed point :
// contains  iCrDistance x iPeriodChild  angles 
// aproximates Julia set near alfa
// target set for forard iterations for all points of interior
double RaysTurns[iMaxDistance2Alfa][iPeriodChild];



static double TwoPi=2*M_PI;

// // http://en.wikipedia.org/wiki/Cohen%E2%80%93Sutherland
// 2d line clipping 
const int INSIDE = 0; // 0000
const int LEFT = 1;   // 0001
const int RIGHT = 2;  // 0010
const int BOTTOM = 4; // 0100
const int TOP = 8;    // 1000

/* ------------------------------------------ functions -------------------------------------------------------------*/

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


// colors of components interior = shades of gray
int InitColors(int iMax, unsigned char iColors[])
{
  int i;
  //int iMax = iPeriodChild; iPeriodChild and iColorsOfInterior
  unsigned int iStep;

  iStep=  2 ;   //150/iMax;

  for (i = 1; i <= iMax; ++i)
    {
      iColors[i-1] = iColorOfExterior -i*iStep; 
      // printf("i= %d color = %i  \n",i-1, iColors[i-1]); // debug
    }
  return 0;
}




//------------------complex numbers -----------------------------------------------------

// gives argument of complex number in turns 
// realted with alfa fixed point
// http://en.wikipedia.org/wiki/Turn_%28geometry%29

double GiveTurn(double zx, double zy)
{
  double argument;

  argument =atan2(zy, zx);// carg(zx +zy*I );   argument in radians from -pi to pi
  if (argument<0) argument=argument + TwoPi; //   argument in radians from 0 to 2*pi
  return argument/TwoPi ; // argument in turns from 0.0 to 1.0
}




/* 
   principal square  root of complex number 
   http://en.wikipedia.org/wiki/Square_root

   z1= I;
   z2 = root(z1);
   printf("zx  = %f \n", creal(z2));
   printf("zy  = %f \n", cimag(z2));
*/
double complex root(double x, double y)
{ 
  
  double u;
  double v;
  double r = sqrt(x*x + y*y); 
  
  v = sqrt(0.5*(r - x));
  if (y < 0) v = -v; 
  u = sqrt(0.5*(r + x));
  return u + v*I;
}



/*
finds offset in turns of parabolic Julia set 

Method is described here :
 Complex Dynamics by Lennart Carleson,Theodore Williams Gamelin 
 page 40

------------Explanation and code by Wolf Jung ------------------------
> The problem is that for q = 10 it is computing 1025 coefficients
> where only 12 are needed. Try the following and compare the
> results with the existing program:
> * Allocate a matrix array a_mn with m <= q and n <= q+1.
> * Set a_11 = lambda , a_12 = 1 , and a_13 ... a_1,q+1 = 0.
> * For m = 1 to q-1 for n = 1 to q+1
>  a_m+1,n = lambda*a_mn + sum a_mk * a_m,n-k
>  where the summation is from k = 1 to n-1.
> Then a_mn is the coefficient of z^n in f^m , with
> f = lambda*z + z^2 . So the desired A is a_q,q+1 .
> 
> In fact you only need one row instead of a quadratic array
> by redefining the variables:
> * Allocate an array a_n with n <= q+1.
> * Set a_1 = lambda , a_2 = 1 , and a_3 ... a_q+1 = 0.
> * For m = 1 to q-1 for n = q+1 down to 1 (not up !)
>  a_n = lambda*a_n + sum a_k * a_n-k
>  where the summation is from k = 1 to n-1.
> 
> Or in c++:
> int k, m, n; complex a[q+2];
> a[1] = lambda; a[2] = complex(1.0);
> for (n = 3; n <= q+1; n++) a[n] = 0;
> for (m = 1; m < q; m++) for (n = q+1; n; n--)
> { a[n] *= lambda;
>  for (k = 1; k < n; k++) a[n] += a[k]*a[n-k];
> }
> //Now a[q+1] = A with f^q(z) = z + a*z^(q+1) + ...
> 
>  Best regards,
> 
>  Wolf

*/
double GiveOffset(int p, int q )
{
  double complex a[iPeriodChild+2];
  double complex lambda;
  double t; // = Internal Angle In radians
  double offset;
  int k, m, n; // index of for loops
  double complex coeff;
  
  t = ((double)p/q) *2*M_PI; // from turns to radians
  lambda= (cos(t)+I*sin(t));
   // a[0]=0;
   a[1] = lambda; 
   a[2] = 1.0;
   for (n = 3; n <= q+1; n++) a[n] = 0;
   //
   for (m = 1; m < q; m++) 
     for (n = q+1; n; n--)
      { a[n] *= lambda;
        for (k = 1; k < n; k++) 
          a[n] += a[k]*a[n-k]; // sum 
      }
   //Now a[q+1] = A with f^q(z) = z + a*z^(q+1) + ...
   coeff = a[q+1];
   offset= -GiveTurn(creal(coeff), cimag(coeff))/q;
  return offset;


}






// from screen to world coordinate ; linear mapping
// uses global cons
double GiveZx(unsigned int ix)
{ return (ZxMin + ix*PixelWidth );}

// uses globaal cons
double GiveZy(unsigned int iy)
{ return (ZyMax - iy*PixelHeight);} // reverse y axis




/* -----------  array functions = drawing -------------- */


/* gives position of 2D point (ix,iy) in 1D array  ; uses also global variable iWidth */
unsigned int Give_i(unsigned int ix, unsigned int iy)
{ return ix + iy*iWidth; }




// plots raster point (ix,iy) 
int iDrawPoint(unsigned char A[], unsigned int ix, unsigned int iy, unsigned char iColor)
{ 

  /* i =  Give_i(ix,iy) compute index of 1D array from indices of 2D array */
  A[Give_i(ix,iy)] = iColor;

  return 0;
}


// draws point to memmory array data
// uses complex type so #include <complex.h> and -lm 
int dDrawPoint(unsigned char A[], complex double point,unsigned char iColor )
{

  unsigned int ix, iy; // screen coordinate = indices of virtual 2D array
  //unsigned int i; // index of 1D array
  
  ix = (creal(point)- ZxMin)/PixelWidth; 
  iy = (ZyMax - cimag(point))/PixelHeight; // inverse Y axis 
  iDrawPoint(A, ix, iy, iColor);
  return 0;
}





// compute turn and save to array
int iSaveTurn(int iX, int iY, int NumberOfRay)
{
  

  double dXt; // = dX-dAlfaX = translation
  double dYt; // = dY-dAlfaY
  double turn; 
  // distance to alfa fixed point ( target set )
  int iDistance;
  double dDistance2; // 
  

  dXt = ( GiveZx(iX) - dAlfaX);
  dYt = ( GiveZy(iY) - dAlfaY );
  dDistance2 = dXt*dXt + dYt*dYt;
  if (dMaxDistance2Alfa2>dDistance2) 
    {

      turn =  GiveTurn(  dXt, dYt);
      iDistance = (int)(sqrt(dDistance2)/PixelWidth);
      // printf("Zx = %f; Zy= %f iDistance = %d\n",GiveZx(iX),GiveZx(iY), iDistance); // debug info 
      if (iDistance<=iMaxDistance2Alfa) RaysTurns[iDistance][NumberOfRay] = turn;
    }
  return 0;

}





/*
  http://rosettacode.org/wiki/Bitmap/Bresenham%27s_line_algorithm
  Instead of swaps in the initialisation use error calculation for both directions x and y simultaneously:
*/
void iDrawLine(unsigned char A[], int x0, int y0, int x1, int y1, int RayNumber, int color) 
{
  int x=x0; int y=y0;
  int dx = abs(x1-x0), sx = x0<x1 ? 1 : -1;
  int dy = abs(y1-y0), sy = y0<y1 ? 1 : -1; 
  int err = (dx>dy ? dx : -dy)/2, e2;
  

  for(;;){
    iDrawPoint(A, x, y,color);
    iSaveTurn( x,y,RayNumber);
    
    if (x==x1 && y==y1) break;
    e2 = err;
    if (e2 >-dx) { err -= dy; x += sx; }
    if (e2 < dy) { err += dx; y += sy; }
  }
}




// compute turn and save to array
// ray goes : from inf thru z0 thru z1 to alfa
// it means that external radius of z0 is greater than external radius of z0
// it does not mean that iDistance0 iDistance1 ( because ray can turn )
double SaveTurns(double Zx0, double Zy0, double Zx1, double Zy1, int NumberOfRay)
{
  double xt; // = x-dAlfaX = translation
  double yt; // = y-dAlfaY
  int iDistance, iDistance0, iDistance1;
  int iDistanceMax; //, iDistanceMin;
  double turn, turn0, turn1; //, turnMin;
  //double dStep; // step between turns
  //int iStep;   // = iDistanceMax - iDistanceMin;

  // compute/save turn of only one point
  // translation to alfa 
  xt= Zx0 - dAlfaX;
  yt= Zy0 - dAlfaY;
     
  iDistance0 = (int)(sqrt(xt*xt +  yt*yt)/PixelWidth);
  turn0 =  GiveTurn(  xt, yt);
  //if inside trap !! save turns to the RaysTurns array
  if ( iDistance0 <= iMaxDistance2Alfa && iDistance0 >0 ) 
    RaysTurns[iDistance0][NumberOfRay]= turn0; 
       

  // compute/save turn of only one point
  // translation to alfa 
  xt= Zx1 - dAlfaX;
  yt= Zy1 - dAlfaY;
      
  iDistance1 = (int)(sqrt(xt*xt +  yt*yt)/PixelWidth); 
  turn1 =  GiveTurn(  xt, yt);
  //if inside trap !! save turns to the RaysTurns array
  if ( iDistance1 <= iMaxDistance2Alfa && iDistance1 >0 ) 
    RaysTurns[iDistance1][NumberOfRay]= turn1; 

       
  // if one of 2 ends is inside trap 
  // then fill the gaps ( where turn <0 ; it means not filled )
  // and save turns to the RaysTurns array
  if ( iDistance0 <= iMaxDistance2Alfa || iDistance1 <= iMaxDistance2Alfa)
    { 
      // all points are not alfa ; both iDistance >0
      //    if ( iDistance0 >0 && iDistance1 > 0  ) 
      /*       {*/
      /*         iStep = iDistanceMax - iDistanceMin;*/
      /*         if (turn0>turn1) */
      /*            {dStep = (turn0-turn1)/iStep; turnMin = turn1;}*/
      /*            else {dStep = (turn1-turn0)/iStep; turnMin=turn0;}*/
      /*         for (iDistance=iMaxDistance2Alfa; iDistance<iDistance1; iDistance++) RaysTurns[iDistance][NumberOfRay] = turnMin+(iDistance-iDistance1)*dStep;*/
      /*         */
      /*       }*/
      /*     */
      // both ends are inside trap 
      // one of 2 points is alfa, so one iDistance ==0
      if (iDistance1==0 ) 
	{iDistanceMax=iDistance0;  turn = turn0;}
      else {iDistanceMax=iDistance1;  turn = turn1;}
      // copy first positive value  to all cells before it = straight line 
      for (iDistance=0; iDistance<iDistanceMax; iDistance++) RaysTurns[iDistance][NumberOfRay] = turn;
        
    }

  return 0;

}

/*

 */
// http://en.wikipedia.org/wiki/Cohen%E2%80%93Sutherland
// Compute the bit code for a point (x, y) using the clip rectangle
// bounded diagonally by (xmin, ymin), and (xmax, ymax)

// ASSUME THAT xmax, xmin, ymax and ymin are global constants.

int ComputeOutCode(double x, double y)
{
  int code;

  code = INSIDE;          // initialised as being inside of clip window

  if (x < ZxMin)  
    code = LEFT; // to the left of clip window
  else if (x > ZxMax) code = RIGHT;     // to the right of clip window
		
  if (y < ZyMin)        // below the clip window
    code = BOTTOM;
  else if (y > ZyMax) code = TOP;     // above the clip window
		

  return code;
}

// http://en.wikipedia.org/wiki/Cohen%E2%80%93Sutherland
// Cohen–Sutherland clipping algorithm clips a line from
// P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with 
// diagonal from (xmin, ymin) to (xmax, ymax).
//void CohenSutherlandLineClipAndDraw(
void dDrawLine(unsigned char A[],double x0, double y0, double x1, double y1, int RayNumber, int color )
{
  // compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
  int outcode0 = ComputeOutCode(x0, y0);
  int outcode1 = ComputeOutCode(x1, y1);
  int accept = 0;
  unsigned int ix0, iy0; // screen coordinate = indices of virtual 2D array 
  unsigned int ix1, iy1; // screen coordinate = indices of virtual 2D array

  while (1) {
    if (!(outcode0 | outcode1)) { // Bitwise OR is 0. Trivially accept and get out of loop
      accept = 1;
      break;
    } else if (outcode0 & outcode1) { // Bitwise AND is not 0. Trivially reject and get out of loop
      break;
    } else {
      // failed both tests, so calculate the line segment to clip
      // from an outside point to an intersection with clip edge
      double x, y;

      // At least one endpoint is outside the clip rectangle; pick it.
      int outcodeOut = outcode0 ? outcode0 : outcode1;

      // Now find the intersection point;
      // use formulas y = y0 + slope * (x - x0), x = x0 + (1 / slope) * (y - y0)
      if (outcodeOut & TOP) {           // point is above the clip rectangle
	x = x0 + (x1 - x0) * (ZyMax - y0) / (y1 - y0);
	y = ZyMax;
      } else if (outcodeOut & BOTTOM) { // point is below the clip rectangle
	x = x0 + (x1 - x0) * (ZyMin - y0) / (y1 - y0);
	y = ZyMin;
      } else if (outcodeOut & RIGHT) {  // point is to the right of clip rectangle
	y = y0 + (y1 - y0) * (ZxMax - x0) / (x1 - x0);
	x = ZxMax;
      } else if (outcodeOut & LEFT) {   // point is to the left of clip rectangle
	y = y0 + (y1 - y0) * (ZxMin - x0) / (x1 - x0);
	x = ZxMin;
      }

      // Now we move outside point to intersection point to clip
      // and get ready for next pass.
      if (outcodeOut == outcode0) {
	x0 = x;
	y0 = y;
	outcode0 = ComputeOutCode(x0, y0);
      } else {
	x1 = x;
	y1 = y;
	outcode1 = ComputeOutCode(x1, y1);
      }
    }
  }



  if (accept) 
    { // draw line clipped to view window 

     
      // SaveTurns(x0,y0,x1,y1,RayNumber);
      // all segment is inside a window 
      ix0= (x0- ZxMin)/PixelWidth; 
      iy0 = (ZyMax - y0)/PixelHeight; // inverse Y axis 
      ix1= (x1- ZxMin)/PixelWidth; 
      iy1= (ZyMax - y1)/PixelHeight; // inverse Y axis 
      iDrawLine(A, ix0,iy0,ix1,iy1,RayNumber,color) ;}
	
}












 




 
int ColourInternalTrap(unsigned char A[])
{
 
  unsigned int ix, iy; // pixel coordinate 
  double Zx, Zy; //  Z= Zx+ZY*i;
  double Zxt, Zyt;
  unsigned i; /* index of 1D array */
  double dDistance2; // = sqr( distance to alfa)
  int iDistance;
  double turn;
   
  for(iy=iyMin;iy<=iyMax;++iy) 
    {
      Zy = GiveZy(iy);
      Zyt = Zy -dAlfaY;// translation near fixed point alfa
      for(ix=ixMin;ix<=ixMax;++ix) 
        {
 
          // from screen to world coordinate 
          Zx = GiveZx(ix);
          
          // translation near fixed point alfa
          Zxt = Zx - dAlfaX; 
          //
          turn =  GiveTurn(Zxt, Zyt);         
          dDistance2 = Zxt*Zxt + Zyt*Zyt ;
          iDistance = (int)(sqrt(dDistance2)/PixelWidth);
          // check if point z is inside triangle around arm of critical orbit 
          if ( dDistance2 < dMaxDistance2Alfa2 && turn>RaysTurns[iDistance][iPeriodChild-2] && turn<RaysTurns[iDistance][iPeriodChild-1]) 
	    { i = Give_i(ix, iy); /* compute index of 1D array from indices of 2D array */
              if ( A[i]==255) A[i]=0; // rays 
	      else A[i]= 255-A[i];   // mark the trap
	    }
        }
    }
 
 
 
  return 0;
}
 
//---------------------------------------------------------------



// save file for "manual" checking it's content'
int SaveFiles4ManualDebug(double A[iMaxDistance2Alfa][iPeriodChild], double dNumberForName)
{

  int d,p;
  char name [30]; /* name of the file */
  sprintf(name,"RaysTurns%f", dNumberForName); /* create name from number */
  char *filename =strcat(name,".txt");
  FILE * fp;

  // save whole A array to the text file
  fp= fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  for (d=0; d<iMaxDistance2Alfa ; ++d)
    { fprintf(fp, " iDistance2Alfa = %d \t", d);
      {  for (p=0; p<iPeriodChild ; ++p) fprintf(fp,"turn(ray%d) = %f \t",p, A[d][p]);
	fprintf(fp, "\n");}  } 
  fclose(fp); 
  //printf(" file %s saved \n", filename);

 
  for (p=0; p<iPeriodChild ; ++p)
    {  // save ray from A array to the text file
      sprintf(name,"Ray%f", p+dNumberForName); /* create name from number */
      filename =strcat(name,".txt");
      fp= fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
      for (d=0; d<iMaxDistance2Alfa ; ++d) fprintf(fp, " [%d ,%f ], ", d,  A[d][p]);
      fclose(fp); 
      //printf(" file %s saved \n", filename);
    }


 


  return 0;


}

/* 
   Maxima CAS code using draw package 
   ( interface to gnuplot ) by Mario Rodriguez Riotorto http://www.telefonica.net/web2/biomates 
   A:[];

*/
/*

  load(draw); 

  draw2d(
  title = "",
  terminal  = screen,
  user_preamble = "set size square", 
  dimensions=[1000,1000],
   
  xlabel     = "z.re ",
  ylabel     = "z.im",
  yrange = [0.0,1.0],
     
  point_type    = filled_circle,
  points_joined = true,
  point_size    = 0.7,
   
  key=" orbit ",
  color             =red,
  points(A)
  );

*/



// aproximate not computed ( negative values)
int FillGapsInRaysArray()
{
     // indexes of the array 
      int i; // index of the pixel on the ray
      int iMini, iMaxi; // gap border
      int j; // index of the ray
 
      double dStep;
      int iStep;

		// fill empty values = change negative with positive !!!!!!!!!!
		// negative value = not filled = gaps
		// positive value = filled ( computed )


		// check  every ray 
		for (j = 0; j < iPeriodChild ; j++)
        {

         // start from 0 wher 
         i=0;
         // rest of the array
         do 
	     {
	      // go thru all positive 
	      while (RaysTurns[i][j]>0.0 && i<iMaxDistance2Alfa-1) i+=1;
	      iMini=i-1; // here is positive, lower border of the gap 
	      // here value is negative : RaysTurns[i][j]<0.0
	      do i+=1; while (RaysTurns[i][j]<0.0 && i<iMaxDistance2Alfa-1  ); 
	      iMaxi= i; // positive , upper border of the gap
	      iStep= iMaxi-iMini;
	      // step between RaysTurns[iMini][j] and RaysTurns[iMaxi][j]
	      if (RaysTurns[iMaxi][j]>RaysTurns[iMini][j] ) 
	          dStep = (RaysTurns[iMaxi][j]-RaysTurns[iMini][j])/iStep;
	          else dStep = (RaysTurns[iMini][j]-RaysTurns[iMaxi][j])/iStep;
	       // step between values 
	      if (iMaxi<iMaxDistance2Alfa-1 && RaysTurns[iMaxi-1][j]< 0.0) // array is numberd from 0 to ...
	        {
	             //printf("in ray j= %d there is a gap from i= %d to i= %d iStep = %d ; dStep = %f\n",j, iMini, iMaxi, iStep, dStep );
	              // fill the gap 
	              i= iMini;
	              while (i<iMaxi-1) { i+=1; RaysTurns[i][j]=RaysTurns[iMini][j] + (i-iMini)*dStep;} 
	              i=iMaxi;
	        }
	       // else {
	       //  printf("in ray j= %d there is a gap from i= %d to i= %d iStep = %d ; \n",j,  iMini, iMaxi, iStep);
	       // copy last positive value  to all cells after it = straight line !!!!
	       //  for (i=iMini+1; i<iMaxi; i++) RaysTurns[i][j] = RaysTurns[iMini][j];
	       //} 
       
         } while (i<iMaxDistance2Alfa);
  
        } // j


       return 0;
}


// This function only works for periodic  angles.
// You must know the iPeriodChild n before calling this function.
// draws all "iPeriodChild" external rays 
// http://commons.wikimedia.org/wiki/File:Backward_Iteration.svg


double complex ComputeRays( unsigned char A[],
			    int n, //iPeriodChild of ray's angle under doubling map
			    int iterMax
			    )
{
  double xNew; // new point of the ray
  double yNew;
  //double xt; // x-dAlfaX
  //double yt; // y-dAlfaY
  double xend ; // re of the endpoint of the ray
  double yend; // im of the endpoint of the ray
  const double R = 10; // very big radius = near infinity
  int j; // number of ray 
  int iter; // index of backward iteration
  double t,t0; // external angle in turns 
  
  //double xt,yt; // xtranslation to alfa
  
  double complex zPrev;
  double u,v; // zPrev = u+v*I
 // double complex zNext;
  //double dDistance;
  //int iDistance ; // dDistance/PixelWidth
  double offset;


  // external angle of the first ray 
  t0 = 1.0/( pow(2.0,iPeriodChild) -1.0); // http://fraktal.republika.pl/mset_external_ray_m.html
  t=t0;


  /* dynamic 1D arrays for coordinates ( x, y) of points with the same R on preperiodic and periodic rays  */
  double *RayXs, *RayYs;
  int iLength = n+2; // length of arrays ?? why +2
  //  creates arrays :  RayXs and RayYs  and checks if it was done
  RayXs = malloc( iLength * sizeof(double) );
  RayYs = malloc( iLength * sizeof(double) );
  if (RayXs == NULL || RayYs==NULL)
    {
      fprintf(stderr,"Could not allocate memory");
      getchar(); 
      return 1; // error
    }
  

//printf("t = %.10f \n", t);
  //  starting points on periodic rays 
  //  with angles t, 2t, 4t...  and the same radius R
  for (j = 0; j < n; j++)
    { // z= R*exp(2*Pi*t)
      RayXs[j] = R*cos((2*M_PI)*t); 
      RayYs[j] = R*sin((2*M_PI)*t);
    
      t *= 2; // t = 2*t
      if (t > 1) t--; // t = t modulo 1 
	  printf("t = %.10f \n", t);
    }
 // zNext = RayXs[0] + RayYs[0] *I;

  // printf("RayXs[0]  = %f \n", RayXs[0]);
  // printf("RayYs[0]  = %f \n", RayYs[0]);

  // z[k] is n-periodic. So it can be defined here explicitly as well.
  RayXs[n] = RayXs[0]; 
  RayYs[n] = RayYs[0];
  

  //   backward iteration of each point z
  for (iter = -10; iter <= iterMax; iter++)
    { 
     	
      for (j = 0; j < n; j++) // iPeriodChild +preperiod
	{ // u+v*i = sqrt(z-c)   backward iteration in fc plane 
	  zPrev = root(RayXs[j+1] - Cx , RayYs[j+1] - Cy ); // , u, v
	  u=creal(zPrev);
	  v=cimag(zPrev);
                
	  // choose one of 2 roots: u+v*i or -u-v*i
	  if (u*RayXs[j] + v*RayYs[j] > 0) 
	    { xNew = u; yNew = v; } // u+v*i
	  else { xNew = -u; yNew = -v; } // -u-v*i
                


	  // draw part of the ray = line from zPrev to zNew
	  dDrawLine(A, RayXs[j], RayYs[j], xNew, yNew, j, 255);
                
	  // check if ray is crossing 0 axis = info for debug ; one can comment it 
	 // if (xNew>dAlfaX && ((yNew>dAlfaY && dAlfaY>RayYs[j]) || (yNew<dAlfaY && dAlfaY<RayYs[j])  ))
	   // {  // translation to alfa 
	  //    xt= xNew - dAlfaX;
	   //   yt= yNew - dAlfaY;
	      //dDistance = sqrt(xt*xt +  yt*yt); 
	      //iDistance = (int)(dDistance/PixelWidth); 
	      //printf("ray %d is crossing 0 axis %f = %d pixels from alfa ; \n",j,dDistance,iDistance); // info for debug
	  //  }
	  //  
	  RayXs[j] = xNew; RayYs[j] = yNew;
                
	       
                  


                
	} // for j ...

          //RayYs[n+k] cannot be constructed as a preimage of RayYs[n+k+1]
      RayXs[n] = RayXs[0]; 
      RayYs[n] = RayYs[0];
          
      // convert to pixel coordinates 
      //  if z  is in window then draw a line from (I,K) to (u,v) = part of ray 
   
      // printf("for iter = %d cabs(z) = %f \n", iter, cabs(RayXs[0] + RayYs[0]*I));
     
    }

	
	
  // Complex Dynamics by Lennart Carleson,Theodore Williams Gamelin  page 40 
  offset=GiveOffset(1,n);
  // aproximate end of ray by straight line to it's landing point here = alfa fixed point
  for (j = 0; j < n + 1; j++)
  {  // compute one new point on the ray using offset
     t = ((double) j)/((double) n) - offset;  // repelling direction ; Complex Dynamics by Lennart Carleson,Theodore Williams Gamelin  page 40 
     xNew = dAlfaX;
     yNew = dAlfaY;
     dDrawLine(A, RayXs[j],RayYs[j], xNew, yNew,j, 255 ); 
	 dDrawLine(A, xNew, yNew, dAlfaX, dAlfaY,j, 255 ); 
  }
 
  // last point of a ray 0
  xend = RayXs[0];
  yend = RayYs[0];

  
  
  // this check can be done only from inside this function
 // t=t0;
 // for (j = 0; j < iPeriodChild + 1; j++)
 //   {
       // iDistance = (int)(sqrt((RayXs[j]-dAlfaX)*(RayXs[j]-dAlfaX) +  (RayYs[j]-dAlfaY)*(RayYs[j]-dAlfaY))/PixelWidth);
        //printf("landing point of ray for angle = %f is = (%f ;  %f ) ; iDistnace = %d \n",t, RayXs[j], RayYs[j],  iDistance);
      //  t *= 2; // t = 2*t
    //} // end of the check 

 


  // free memmory
  free(RayXs);
  free(RayYs);





  SaveFiles4ManualDebug(RaysTurns, 0.1); // save filled array with gaps ( turn = -1.000)
  FillGapsInRaysArray();
  // ComputeMinPosition(); // only after FillGapsInRaysArray
  //SaveFiles4ManualDebug(RaysTurns, 0.2); // save filled array without gaps ( all 0<= turns < 1)

  

  return  xend + yend*I; // return last point or ray for angle t 
}







//;;;;;;;;;;;;;;;;;;;;;;  setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

int setup(int ParentPeriod, int ChildPeriod)
{

  int iDistance;
  int j; 
  unsigned int denominator;
  double InternalAngle;
  
  printf("setup\n");

  denominator = ChildPeriod;
  InternalAngle = ((double) numerator)/((double) denominator);

  c = GiveC(InternalAngle, 1.0, ParentPeriod) ; // internal radius= 1.0 gives root point = parabolic parameter 
  Cx=creal(c);
  Cy=cimag(c);
  alfa = GiveAlfaFixedPoint(c);
  dAlfaX = creal(alfa);
  dAlfaY = cimag(alfa);
  iAlfaX = (dAlfaX- ZxMin)/PixelWidth; 
  iAlfaY = (ZyMax - dAlfaY)/PixelHeight; // inverse Y axis 
  
  

  /* 2D array ranges */
  if (!(iHeight % 2)) iHeight+=1; // it sholud be even number (variable % 2) or (variable & 1)
  iWidth = iHeight;
  iSize = iWidth*iHeight; // size = number of points in array 
  // iy
  iyMax = iHeight - 1 ; // Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].
  //ix
  
  ixMax = iWidth - 1;

  /* 1D array ranges */
  // i1Dsize = i2Dsize; // 1D array with the same size as 2D array
  iMax = iSize-1; // Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].


  /* Pixel sizes */
  PixelWidth = (ZxMax-ZxMin)/ixMax; //  ixMax = (iWidth-1)  step between pixels in world coordinate 
  PixelHeight = (ZyMax-ZyMin)/iyMax;
  ratio = ((ZxMax-ZxMin)/(ZyMax-ZyMin))/((float)iWidth/(float)iHeight); // it should be 1.000 ...
  
  

  // for numerical optimisation in iteration
  ER2 = ER * ER;
  iMaxDistance2Alfa2 =iMaxDistance2Alfa * iMaxDistance2Alfa;
  dMaxDistance2Alfa2 = iMaxDistance2Alfa2*PixelWidth*PixelWidth; // dMaxDistance2Alfa^2
  dMaxDistance2Alfa = sqrt(dMaxDistance2Alfa2); // maybe it should be in reversed order ??
    
  
  /* create dynamic 1D arrays for colors ( shades of gray ) */
  data = malloc( iSize * sizeof(unsigned char) );
  rays = malloc( iSize * sizeof(unsigned char) );
  edge = malloc( iSize * sizeof(unsigned char) );

  if (rays == NULL || edge== NULL || data == NULL)
    {
      fprintf(stderr," Could not allocate memory");
      getchar(); 
      return 1;
    }

  // fill array g with negative number 
  for (iDistance=0; iDistance<iMaxDistance2Alfa; ++iDistance)
    for (j=0; j<iPeriodChild; ++j)
      RaysTurns[iDistance][j]=-1.0;

  // fill array iColorsOfInterior with iPeriodChild colors ( shades of gray )
  InitColors(iPeriodChild, iColorsOfInterior);


 
  

  
   
  
  printf(" end of setup \n");
  
  return 0;

} // ;;;;;;;;;;;;;;;;;;;;;;;;; end of the setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



unsigned char ComputeColor(unsigned int ix, unsigned int iy, int IterationMax)
   { 
     // check behavour of z under fc(z)=z^2+c
     // using 2 target set:
     // 1. exterior or circle (center at origin and radius ER ) 
     // as a target set containing infinity = for escaping points ( bailout test)
     // for points of exterior of julia set
	 // 2. interior of circle with center = alfa and radius dMaxDistance2Alfa
	 // as a target set for points of interior of Julia set 
     //  Z= Zx+ZY*i;

     double Zx2, Zy2;
     int i=0;
     //int j; // iteration = fc(z)
     double d2 ; /* d2= (distance from z to Alpha)^2   */
     double Zxt,Zyt ; // 
     double Zx, Zy;
     double turn;
     int iDistance;
  
  
     // from screen to world coordinate 
     Zx = GiveZx(ix);
     Zy = GiveZy(iy);
     /* distance from z to Alpha  */
     Zxt=Zx-dAlfaX;
     Zyt=Zy-dAlfaY;
     d2=Zxt*Zxt +Zyt*Zyt;
     if (d2<dMaxDistance2Alfa2) 
        {
           iDistance = (int)(sqrt(d2)/PixelWidth);
           if (iDistance<iMaxDistance2Alfa)
	        { 
	          turn =  GiveTurn(Zxt, Zyt);
	          if ( turn>RaysTurns[iDistance][iPeriodChild-2] && turn<RaysTurns[iDistance][iPeriodChild-1])
	          return  iColorsOfInterior[i % iPeriodChild];
	        } 
        }

      // if not inside target set around attractor ( alfa fixed point )
       while (1 )
        { // then iterate 
      
          Zx2 = Zx*Zx; 
          Zy2 = Zy*Zy;
       
      // bailout test 
      if (Zx2 + Zy2 > ER2) return iColorOfExterior; // if escaping stop iteration
       
      // if not escaping or not attracting then iterate = check behaviour
      // new z : Z(n+1) = Zn * Zn  + C
      Zy = 2*Zx*Zy + Cy; 
      Zx = Zx2 - Zy2 + Cx; 
      //
      i+=1;
	 
      /* distance from z to Alpha  */
      Zxt=Zx-dAlfaX;
      Zyt=Zy-dAlfaY;
      d2=Zxt*Zxt +Zyt*Zyt;
      // check if fall into internal target set 
      if (d2<dMaxDistance2Alfa2) 
	   {
	      iDistance = (int)(sqrt(d2)/PixelWidth);
	      if (iDistance<iMaxDistance2Alfa)
	       { 
	         turn =  GiveTurn(Zxt, Zyt);
	         if ( turn>RaysTurns[iDistance][iPeriodChild-2] && turn<RaysTurns[iDistance][iPeriodChild-1])
		      return  iColorsOfInterior[i % iPeriodChild];
	        }
	    }
    
      if (i > IterationMax) return 255;
      
      
    }

  return  iColorsOfInterior[i % iPeriodChild];   //
}




// plots raster point (ix,iy) 
int PlotPoint(unsigned char A[] , unsigned int ix, unsigned int iy, int IterationMax)
{
  unsigned i; /* index of 1D array */
  unsigned char iColor;
  


  i = Give_i(ix,iy); /* compute index of 1D array from indices of 2D array */
  iColor = ComputeColor(ix, iy, IterationMax);
  A[i] = iColor;

  return 0;
}


// fill array 
// uses global var :  ...
// scanning complex plane 
int ComputeFatouComponents(unsigned char A[], int IterationMax )
{
  unsigned int ix, iy; // pixel coordinate 

  //printf("compute image \n");
  // for all pixels of image 
  for(iy = iyMin; iy<=iyMax; ++iy) 
    { //printf(" %d z %d\n", iy, iyMax); //info 
      for(ix= ixMin; ix<=ixMax; ++ix) PlotPoint(A, ix, iy, IterationMax ) ; //  
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


int CopyRaysTo(unsigned char A[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
 
 
  printf("copy boundaries from edge array to data array \n");
  for(iY=1;iY<iyMax-1;++iY)
    for(iX=1;iX<ixMax-1;++iX)
      {i= Give_i(iX,iY); if (rays[i]==255) A[i]=0;}
 
 
 
  return 0;
}


 
 
int DrawCriticalOrbit(unsigned char A[], unsigned int IterMax)
{
  // integer = pixel coordinate
  unsigned int ix, iy; 
  // double = world coordinate
  // critical point Z= Zx+ZY*i;
  double Zx=0.0;
  double Zy=0.0; 
  double Zx2=0.0;
  double Zy2=0.0;
  
  unsigned int i; /* index of 1D array */
  unsigned int j; // number of iteration
 
 
  // draw critical point  
  ix = (int)((Zx-ZxMin)/PixelWidth);
  iy = (int)((ZyMax-Zy)/PixelHeight); // reverse y axis
  i = Give_i(ix, iy); /* compute index of 1D array from indices of 2D array */
  A[i]=255-A[i]; 
 
  // iterate
  for (j = 1; j <= IterMax; j++) //larg number of iteration s
    {  Zx2 = Zx*Zx; 
      Zy2 = Zy*Zy;
 
      // bailout test 
      if (Zx2 + Zy2 > ER2) return iColorOfExterior; // if escaping stop iteration
 
      // if not escaping iterate
      // Z(n+1) = Zn * Zn  + C
      Zy = 2*Zx*Zy + Cy; 
      Zx = Zx2 - Zy2 + Cx;
      //compute integer coordinate  
      ix = (int)((Zx-ZxMin)/PixelWidth);
      iy = (int)((ZyMax-Zy)/PixelHeight); // reverse y axis
      i = Give_i(ix, iy); /* compute index of 1D array from indices of 2D array */
      A[i]=255-A[i];   // mark the critical orbit
 
    }
  return 0;
}


// Check Orientation of image : mark first quadrant 
// it should be in the upper right position
// uses global var :  ...
int CheckOrientation(unsigned char A[] )
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
	  if (Zx>0 && Zy>0) A[i]=255-A[i];   // check the orientation of Z-plane by marking first quadrant */

	}
    }
   
  return 0;
}

 


// save "A" array to pgm file 
int SaveArray2PGMFile( unsigned char A[], double k)
{
  
  FILE * fp;
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  char name [30]; /* name of file */
  sprintf(name,"%f", k); /*  */
  char *filename =strcat(name,".pgm");
  char *comment="# ";/* comment should start with # */

  /* save image to the pgm file  */      
  fp= fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P5\n %s\n %u %u\n %u\n",comment,iWidth,iHeight,MaxColorComponentValue);  /*write header to the file*/
  fwrite(A,iSize,1,fp);  /*write A array to the file in one step */
  printf("File %s saved. \n", filename);
  fclose(fp);

  return 0;
}


int info()
{
  // diplay info messages
  printf("Numerical approximation of parabolic Julia set for fc(z)= z^2 + c \n");
  printf("iPeriodParent = %d \n", iPeriodParent);
  printf("iPeriodOfChild  = %d \n", iPeriodChild);
  printf("parameter c = ( %f ; %f ) \n", Cx, Cy); 
  printf("alfa fixed point z = ( %f ; %f )  \n", dAlfaX, dAlfaY);
  printf("Image Width = %f \n", ZxMax-ZxMin);
  printf("PixelWidth = %f \n", PixelWidth);
  printf("size of target set in screen units = iMaxDistance2Alfa  = %d pixels \n", iMaxDistance2Alfa); 
  printf("size of target set in world units = dMaxDistance2Alfa  = %f ; \n", dMaxDistance2Alfa);
  printf("Maximal number of iterations = iterMax = %ld \n", iterMax);
  printf("ratio of image  = %f ; it should be 1.000 ...\n", ratio);
  return 0;
}



/* -----------------------------------------  main   -------------------------------------------------------------*/
int main()
{
  setup(iPeriodParent, iPeriodChild);
 

 
   
  printf(" compute and draw external rays landing on the fixed point alfa \n");
  ComputeRays( rays, iPeriodChild, 10*iterMax);
  SaveArray2PGMFile( rays, iPeriodChild*1000+iMaxDistance2Alfa+0.0); // only rays  to pgm file 
 

  ComputeFatouComponents(data, iterMax);
  SaveArray2PGMFile( data, iPeriodChild*1000+iMaxDistance2Alfa+0.1); // save array data (components of Fatou set ) to pgm file


  ComputeBoundariesIn(data);
  SaveArray2PGMFile( edge, iPeriodChild*1000+iMaxDistance2Alfa+0.2); // save array edge (Julia set ) to pgm file


  CopyRaysTo(data);
  SaveArray2PGMFile( data, iPeriodChild*1000+iMaxDistance2Alfa+0.3); // save array data (components + rays) to pgm file

  CopyRaysTo(edge);
  SaveArray2PGMFile( edge, iPeriodChild*1000+iMaxDistance2Alfa+0.4); // save array (Julia set + rays) to pgm file

  CopyBoundariesTo(data);
  SaveArray2PGMFile( data, iPeriodChild*1000+iMaxDistance2Alfa+0.5); // save array data (Julia set and components ) to pgm file

  

  
  // 
  ColourInternalTrap(rays);
  SaveArray2PGMFile( rays, iPeriodChild*1000+iMaxDistance2Alfa+0.6); // save array rays ( = rays + axis + trap) to pgm file
  

  

  ColourInternalTrap(data);
  SaveArray2PGMFile( data, iPeriodChild*1000+iMaxDistance2Alfa+0.7); // save array data (image  = + trap ) to pgm file

  DrawCriticalOrbit(data, 100000);
  SaveArray2PGMFile(data, iPeriodChild*1000+iMaxDistance2Alfa+0.8); // save array rays ( = rays + axis + trap + critical orbit) to pgm file

  

  printf(" always free memory  to avoid buffer overflow \n");
  free(rays);
  free(data);
  free(edge);

  
  info();

  return 0;
}
