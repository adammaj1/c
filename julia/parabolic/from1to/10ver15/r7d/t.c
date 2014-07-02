/*

  Adam Majewski
  fraktal.republika.pl

  c console progam using 
  * symmetry
  * openMP

  draw  parabolic Julia set 

  and saves it to pgm files ( different versions )

 

  gcc t.c -lm -Wall -fopenmp -march=native 
  time ./a.out


  main parameters :
  - iPeriodChild ( c is a root point between period 1 and period = iPeriodChild components of Mandelbrot set) It means that c is
  from elephant valley

  - FillRaysArray(10000000); number here show how smooth will be a boundary ( Julia set ) near parabolic fixed points
  Of course more smooth means more time to compute iter

  -  iMaxDistance2Alfa = radius of circle around alfa. This circle is a target set ( trap) for points that fall 
  into alfa fixed point  = points of interior of Filled Julia set ( more radius = faster but also maybe distorted image, check it )




*/



#include <stdio.h>
#include <stdlib.h> // malloc
#include <string.h> // strcat
#include <math.h> // M_PI; needs -lm also 
#include <complex.h>
#include <omp.h> // OpenMP; needs also -fopenmp


/* --------------------------------- global variables and constans ------------------------------------------------------------ */
// iPeriodChild of secondary component joined by root point
#define iPeriodChild 15

//  unsigned int denominator;  denominator = iPeriodChild;
double InternalAngle;
unsigned char Colors[iPeriodChild]; //={255,230,180, 160,140,120,100}; // NumberOfPetal of colors = iPeriodChild
unsigned char iExterior = 245;
unsigned char iPetal = 255;

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
unsigned int iHeight = 2000; //  odd NumberOfPetal !!!!!! = (iyMax -iyMin + 1) = iyAboveAxisLength + iyBelowAxisLength +1
// The size of array has to be a positive constant integer 
unsigned int iSize ; // = iWidth*iHeight; 


// memmory 1D arrays 
unsigned char *data;
unsigned char *edge;


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
 
double ratio ;




// complex numbers of parametr plane 
double Cx; // c =Cx +Cy * i
double Cy;
double complex c; // 

double complex alfa; // alfa fixed point alfa=f(alfa)
double alfax,alfay;

unsigned long int iterMax  =100000000; //iHeight*100;
// target set for escaping points is a exterior of circle with center in origin
double ER = 2.0; // Escape Radius for bailout test 
double ER2;



#define iMaxDistance2Alfa  120 // distance point to alfa fixed point in pixels where PixelWidth = 0.003; // 3 world units / 1000 pixels 
double dMaxDistance2Alfa2; // = (iMaxDistance2Alfa*PixelWidth)^2

//
// target set for points falling into alfa fixed point
// is a circle around alfa fixed point
// with radius = AR
//double AR ; // radius of target set around alfa fixed point in world coordinate AR = PixelWidth*TargetWidth; 
//double AR2; // =AR*AR;


double TwoPi=2*M_PI;

// array with angles in turns of points of periodic rays landing on alfa fixed point :
// contains  iCrDistance x iPeriodChild  angles 
double RaysTurns[iMaxDistance2Alfa][iPeriodChild];









/* ------------------------------------------ functions -------------------------------------------------------------*/





// gives argument of complex number in turns 
// realted with alfa fixed point
// http://en.wikipedia.org/wiki/Turn_%28geometry%29
double GiveTurn(double zx, double zy)
{
  double argument;

  argument =atan2(zy-alfay, zx-alfax);// carg(zx-alfax +(zy-alfay)*I ); // alfa !!!;   argument in radians from -pi to pi
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


// This function  works for periodic angles.
// You must determine the period n before calling this function.
// draws all "period" external rays 

// based on the code for backward iteration for drawing external ray see QmnPlane::backRay()
// by Wolf Jung  http://www.mndynamics.com/indexp.html


double complex FillRaysArray(int IterMax  )
{
  
  double xend ; // re of the endpoint of the ray
  double yend; // im of the endpoint of the ray
  const double R = 10000; // very big radius = near infinity
  int j; // number of ray 
  int iter; // index of backward iteration
  double t,t0;
  
  
  double complex zPrev;
  double u,v; // zPrev = u+v*I
 //double complex zNext;

  int iDistance ; // dDistance/PixelWidth
  
  // fill array with negative number 
  for (iDistance=0; iDistance<iMaxDistance2Alfa; ++iDistance)
    for (j=0; j<iPeriodChild; ++j)
      RaysTurns[iDistance][j]=-1.0;


  t0 = 1.0/( pow(2.0,iPeriodChild) -1.0); // http://fraktal.republika.pl/mset_external_ray_m.html
  t=t0;

  /* dynamic 1D arrays for coordinates ( x, y) of points with the same R on preperiodic and periodic rays  */
  double *RayXs, *RayYs;
  int iLength = iPeriodChild+2; // length of arrays ?? why +2
  //  creates arrays :  RayXs and RayYs  and checks if it was done
  RayXs = malloc( iLength * sizeof(double) );
  RayYs = malloc( iLength * sizeof(double) );
  if (RayXs == NULL || RayYs==NULL)
    {
      fprintf(stderr,"Could not allocate memory");
      getchar(); 
      return 1; // error
    }
  


  //  starting points on preperiodic and periodic rays 
  //  with angles t, 2t, 4t...  and the same radius R
  for (j = 0; j < iPeriodChild ; j++)
    { // z= R*exp(2*Pi*t)
      RayXs[j] = R*cos((2*M_PI)*t); 
      RayYs[j] = R*sin((2*M_PI)*t);


    
      t *= 2; // t = 2*t
      if (t > 1) t--; // t = t modulo 1 
    }
 // zNext = RayXs[0] + RayYs[0] *I;

  // ???
  // z[k] is n-periodic. So it can be defined here explicitly as well.
  RayXs[iPeriodChild] = RayXs[0]; 
  RayYs[iPeriodChild] = RayYs[0];
  

  //   backward iteration of each point z
  for (iter = -10; iter <= IterMax; iter++)
    { 
     	
      for (j = 0; j < iPeriodChild; j++) // period +preperiod
	{ // u+v*i = sqrt(z-c)   backward iteration in fc plane 
	  zPrev = root(RayXs[j+1] - creal(c)+(RayYs[j+1] - cimag(c))*I ); // , u, v
	  u=creal(zPrev);
	  v=cimag(zPrev);
                
	  // choose one of 2 roots: u+v*i or -u-v*i
	  if (u*RayXs[j] + v*RayYs[j] > 0) 
	    { RayXs[j] = u; RayYs[j] = v; } // u+v*i
	  else { RayXs[j] = -u; RayYs[j] = -v; } // -u-v*i


	  //if inside trap !! save turns to the array
	  iDistance = (int)(sqrt((RayXs[j]-alfax)*(RayXs[j]-alfax) +  (RayYs[j]-alfay)*(RayYs[j]-alfay))/PixelWidth);
	  if ( iDistance < iMaxDistance2Alfa )
                   
	    {
                   
	      RaysTurns[iDistance][j]= GiveTurn(  RayXs[j],  RayYs[j]);
	    }
                  
                 
                
	} // for j ...


      // ???
      // z[k] is n-periodic. So it can be defined here explicitly as well.
      RayXs[iPeriodChild] = RayXs[0]; 
      RayYs[iPeriodChild] = RayYs[0];   
          
     
    }



  // check 
  t=t0;
  for (j = 0; j < iPeriodChild + 1; j++)
    {
      // aproximate end of ray by straight line to it's landing point here = alfa fixed point
      //dDrawLine(RayXs[j],RayYs[j], creal(alfa), cimag(alfa), 0, data); 
      iDistance = (int)(sqrt((RayXs[j]-alfax)*(RayXs[j]-alfax) +  (RayYs[j]-alfay)*(RayYs[j]-alfay))/PixelWidth);
      printf("landing point of ray for angle = %f is = (%f ;  %f ) ; iDistnace = %d \n",t, RayXs[j], RayYs[j],  iDistance);
      t *= 2; // t = 2*t
    } // end of the check 


  // last point of a ray 0
  xend = RayXs[0];
  yend = RayYs[0];

  

  // free memmory
  free(RayXs);
  free(RayYs);

  return  xend + yend*I; // return last point or ray for angle t 
}







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

      // start from 0 wher are negative values 
      // because of slow dynamics 
      i=0;
      // find first positive value 
      while (RaysTurns[i][j]<0.0) i+=1;
      iMaxi = i;// here positive

      // copy first positive value  to all cells before it = straight line 
      for (i=0; i<iMaxi; i++) RaysTurns[i][j] = RaysTurns[iMaxi][j];
      // go back 
      i = iMaxi; // positive
      printf("in ray j= %d there is a gap from i= 0 to i= %d iStep = %d ; \n",j,  iMaxi, iMaxi);
      // rest of the array
      do {
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
	    printf("in ray j= %d there is a gap from i= %d to i= %d iStep = %d ; dStep = %f\n",j, iMini, iMaxi, iStep, dStep );
	    // fill the gap 
	    i= iMini;
	    // while (i<iMaxi-1) { i+=1; RaysTurns[i][j]=RaysTurns[iMini][j] + (i-iMini)*dStep;} 
	    i=iMaxi;
	  }
        else {
	  printf("in ray j= %d there is a gap from i= %d to i= %d iStep = %d ; \n",j,  iMini, iMaxi, iStep);
	  // copy last positive value  to all cells after it = straight line 
	  for (i=iMini+1; i<iMaxi; i++) RaysTurns[i][j] = RaysTurns[iMini][j];
	} 
       
      } while (i<iMaxDistance2Alfa);
  
    } // j


  return 0;
}


// row of array RaysTurns [j] 
// will be ordered from :
// lowest angle  RaysTurns[j][0]>0.0
// to maximal angle RaysTurns[j][iMaxDistance2Alfa-1] < 1.0
int OrderRaysArray()

{ 

  double row[iPeriodChild];
  double SmallestAngle;
  int i,j;
  
  int iSmallest;

  // find smallest angle in Turns array 
  SmallestAngle = RaysTurns[0][0];
  for (i=0; i<iPeriodChild; ++i) if (RaysTurns[0][i]<SmallestAngle) {SmallestAngle=RaysTurns[0][i]; iSmallest=i;}


  for (i=0; i<iMaxDistance2Alfa; ++i)  

    {
      for (j = 0; j < iPeriodChild ; j++)
	row[j]= RaysTurns[i][(iSmallest+j) % iPeriodChild]; // copy to row arrays
      for (j = 0; j < iPeriodChild ; j++)
	// copy from row to RayTurns   
	RaysTurns[i][j] = row[j];
   

    }
 



  return 0;
}

int SaveArray2File(double a[iMaxDistance2Alfa][iPeriodChild])
{

  int d,p;


  FILE * fp;
  fp= fopen("file.txt","wb"); /*create new file,give it a name and open it in binary mode  */
  for (d=0; d<iMaxDistance2Alfa ; ++d)
    { fprintf(fp, " i= %d ", d);
      {  for (p=0; p<iPeriodChild ; ++p) fprintf(fp,"turn = %f ", a[d][p]);
	fprintf(fp, "\n");}  } 
  fclose(fp); 


  return 0;


}



/* find c in component of Mandelbrot set 
 
   uses code by Wolf Jung from program Mandel
   see function mndlbrot::bifurcate from mandelbrot.cpp
   http://www.mndynamics.com/indexp.html

*/
double complex GiveC(double InternalAngleInTurns, double InternalRadius, unsigned int iPeriod)
{
  //0 <= InternalRay<= 1
  //0 <= InternalAngleInTurns <=1
  double t = InternalAngleInTurns *2*M_PI; // from turns to radians
  double R2 = InternalRadius * InternalRadius;
  //double Cx, Cy; /* C = Cx+Cy*i */
  switch ( iPeriod ) // of component 
    {
    case 1: // main cardioid
      Cx = (cos(t)*InternalRadius)/2-(cos(2*t)*R2)/4; 
      Cy = (sin(t)*InternalRadius)/2-(sin(2*t)*R2)/4; 
      break;
    case 2: // only one component 
      Cx = InternalRadius * 0.25*cos(t) - 1.0;
      Cy = InternalRadius * 0.25*sin(t); 
      break;
      // for each iPeriodChild  there are 2^(iPeriodChild-1) roots. 
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
  ax^2 +bx + c =0 // general form  of quadratic equation
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
int InitColors()
{
  int i;
  int iMax = iPeriodChild; // uses global var iPeriodChild and Colors
  unsigned int iStep;

  iStep=150/iPeriodChild;

  for (i = 1; i <= iMax; ++i)
    {Colors[i-1] = iExterior -i*iStep; 
      printf("i= %d color = %i  \n",i-1, Colors[i-1]);}
  return 0;
}


/* -------------------------------------------------- SETUP --------------------------------------- */
int setup()

{

  
  /* 2D array ranges */
  if (!(iHeight % 2)) iHeight+=1; // it sholud be even NumberOfPetal (variable % 2) or (variable & 1)
  iWidth = iHeight;
  iSize = iWidth*iHeight; // size = NumberOfPetal of points in array 
  
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
   
  // for numerical optimisation 
  ER2 = ER * ER;
  dMaxDistance2Alfa2 = (iMaxDistance2Alfa*PixelWidth)*(iMaxDistance2Alfa*PixelWidth);// AR2 
   
  /* create dynamic 1D arrays for colors ( shades of gray ) */
  data = malloc( iSize * sizeof(unsigned char) );
  edge = malloc( iSize * sizeof(unsigned char) );
  if (edge == NULL || edge == NULL )
    {
      fprintf(stderr," Could not allocate memory\n");
      return 1;
    }
  else fprintf(stderr," memory is OK \n");

  
  InitColors();
  

  // 
  InternalAngle = 1.0/((double) iPeriodChild); // 
  c = GiveC(InternalAngle, 1.0, 1) ;
  Cx=creal(c);
  Cy=cimag(c);
  //
  alfa = GiveAlfaFixedPoint(c);
  alfax=creal(alfa);
  alfay=cimag(alfa);

  // array of turns 
  
  FillRaysArray(iterMax);
  //printf(" dist2 alfa = %d \n", (int)(((creal(z)-alfax)*(creal(z)-alfax) +  (cimag(z)-alfay)*(cimag(z)-alfay))/PixelWidth));
  FillGapsInRaysArray(); 
  OrderRaysArray();
  SaveArray2File(RaysTurns);
  
 
  return 0;

}






// from screen to world coordinate ; linear mapping
// uses global cons
double GiveZx(unsigned int ix)
{ return (ZxMin + ix*PixelWidth );}

// uses global cons
double GiveZy(unsigned int iy)
{ return (ZyMax - iy*PixelHeight);} // reverse y axis










// all points of interior fall into parabolic fixed point z=alfa
// thru iPeriodChild petals
// boundaries of petals are aproximated by periodic rays
// landing on the alfa fixed point 
unsigned char GiveColorOfInterior(double x, double y, double distance2alfa2)
 
{
  double angle;
  int iDistance2Alfa;

 
  iDistance2Alfa = (int)(sqrt(distance2alfa2)/PixelWidth); 
 
  int i;
  

  
  //  
  angle=GiveTurn(x,y); // z-alfa !!!!
  if (angle<RaysTurns[iDistance2Alfa][0] || angle>RaysTurns[iDistance2Alfa][iPeriodChild-1])
    return Colors[0]; 
  for(i=1;i<iPeriodChild-1;++i)  if (angle<RaysTurns[iDistance2Alfa][i]) return Colors[i];

  /*j=0;*/
  /*  //printf("i= %d\n",i);*/
  /*   while (angle > RaysTurns[iDistance2Alfa][j] ) */
  /*    j+=1; */
  return Colors[iPeriodChild-1];  
 
    
 
}



unsigned char GiveColor(unsigned int ix, unsigned int iy)
{ 
  // check behavour of z under fc(z)=z^2+c
  // using 2 target set:
  // 1. exterior or circle (center at origin and radius ER ) 
  // as a target set containing infinity = for escaping points ( bailout test)
  // for points of exterior of julia set
  // 2. interior of circle with center = alfa and radius dMaxDistance2Alfa
  // as a target set for points of interior of Julia set 
  double Zx, Zy; //  Z= Zx+ZY*i;
  double Zx2, Zy2;
  int i;
  //int j; // iteration = fc(z)
  double d2 ; /* d2= (distance from z to Alpha)^2   */
  double dX,dY ; // d = dx+dy*I
  
  // from screen to world coordinate 
  Zx = GiveZx(ix);
  Zy = GiveZy(iy);
  /* distance from z to Alpha  */
  dX=Zx-alfax;
  dY=Zy-alfay;
  d2=dX*dX+dY*dY;
  // if inside target set around attractor ( alfa fixed point )
  while (d2>dMaxDistance2Alfa2)
    {
      for(i=0;i<iPeriodChild ;++i) // iMax = period !!!!
	{ 
	  Zx2 = Zx*Zx; 
	  Zy2 = Zy*Zy;
       
	  // bailout test 
	  if (Zx2 + Zy2 > ER2) return iExterior; // if escaping stop iteration
       
	  // if not escaping or not attracting then iterate = check behaviour
	  // new z : Z(n+1) = Zn * Zn  + C
	  Zy = 2*Zx*Zy + Cy; 
	  Zx = Zx2 - Zy2 + Cx; 
	} 
      /* distance from z to Alpha  */
      dX=Zx-alfax;
      dY=Zy-alfay;
      d2=dX*dX+dY*dY;
      // if inside target set around attractor ( alfa fixed point )
      
      
    }
  return GiveColorOfInterior( Zx, Zy, d2);//iPetal; // not escaping and not in attracting target set  , probably never here (:-))
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
  for(ix=ixMin;ix<=ixMax;++ix) PlotPoint(ix, iy, GiveColor(ix, iy));


  /*
    The use of ‘shared(variable, variable2) specifies that these variables should be shared among all the threads.
    The use of ‘private(variable, variable2)’ specifies that these variables should have a seperate instance in each thread.
  */

#pragma omp parallel for schedule(dynamic) private(iyAbove,ix,iy,Color) shared(iyAboveMin, iyAboveMax,ixMin,ixMax, iyAxisOfSymmetry)

  // above and below axis 
  for(iyAbove = iyAboveMin; iyAbove<=iyAboveMax; ++iyAbove) 
    {// printf(" %d from %d\r", iyAbove, iyAboveMax); //info 
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

int AddBoundaries(unsigned char data[])
{

  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
  /* sobel filter */
  unsigned char G, Gh, Gv; 
 
  


  printf(" find boundaries in data array using  Sobel filter\n");   
#pragma omp parallel for schedule(dynamic) private(i,iY,iX,Gv,Gh,G) shared(iyMax,ixMax, ER2)
  for(iY=1;iY<iyMax-1;++iY){ 
    for(iX=1;iX<ixMax-1;++iX){ 
      Gv= data[Give_i(iX-1,iY+1)] + 2*data[Give_i(iX,iY+1)] + data[Give_i(iX-1,iY+1)] - data[Give_i(iX-1,iY-1)] - 2*data[Give_i(iX-1,iY)] - data[Give_i(iX+1,iY-1)];
      Gh= data[Give_i(iX+1,iY+1)] + 2*data[Give_i(iX+1,iY)] + data[Give_i(iX-1,iY-1)] - data[Give_i(iX+1,iY-1)] - 2*data[Give_i(iX-1,iY)] - data[Give_i(iX-1,iY-1)];
      G = sqrt(Gh*Gh + Gv*Gv);
      i= Give_i(iX,iY); /* compute index of 1D array from indices of 2D array */
      if (G==0) {edge[i]=255;} /* background */
      else {edge[i]=0;}  /* boundary */
    }
  }

  // copy boundaries from edge array to data array 
  //for(iY=1;iY<iyMax-1;++iY){ 
  //   for(iX=1;iX<ixMax-1;++iX){i= Give_i(iX,iY); if (edge[i]==0) data[i]=0;}}



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

int CopyBoundaries()
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
 
 
  printf("copy boundaries from edge array to data array \n");
  for(iY=1;iY<iyMax-1;++iY)
    for(iX=1;iX<ixMax-1;++iX)
      {i= Give_i(iX,iY); if (edge[i]==0) data[i]=0;}
 
 
 
  return 0;
}




int DrawCriticalOrbit(unsigned int IterMax)
{
 
  unsigned int ix, iy; // pixel coordinate 
  double Zx=0.0;
  double Zy=0.0; //  Z= Zx+ZY*i;
  double Zx2=0.0;
  double Zy2=0.0;
  unsigned int i; /* index of 1D array */
  unsigned int j;


  // draw critical point  
  ix = (int)((Zx-ZxMin)/PixelWidth);
  iy = (int)((ZyMax-Zy)/PixelHeight); // reverse y axis
  i = Give_i(ix, iy); /* compute index of 1D array from indices of 2D array */
  data[i]=255-data[i]; 

  // iterate
  for (j = 1; j <= IterMax; j++) //larg number of iteration s
    {  Zx2 = Zx*Zx; 
      Zy2 = Zy*Zy;
       
      // bailout test 
      if (Zx2 + Zy2 > ER2) return iExterior; // if escaping stop iteration
       
      // if not escaping iterate
      // Z(n+1) = Zn * Zn  + C
      Zy = 2*Zx*Zy + Cy; 
      Zx = Zx2 - Zy2 + Cx;
      //compute integer coordinate  
      ix = (int)((Zx-ZxMin)/PixelWidth);
      iy = (int)((ZyMax-Zy)/PixelHeight); // reverse y axis
      i = Give_i(ix, iy); /* compute index of 1D array from indices of 2D array */
      data[i]=255-data[i];   // mark the trap

    }
  return 0;
}





// save data array to pgm file 
int SaveArray2PGMFile( unsigned char data[], double t)
{
  
  FILE * fp;
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  char name [10]; /* name of file */
  sprintf(name,"%f", t); /*  */
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
  printf("InternalAngle  = %f \n", InternalAngle);
  printf("Cx  = %f \n", Cx); 
  printf("Cy  = %f \n", Cy);
  // 
  printf("alfax  = %f \n", creal(alfa));
  printf("alfay  = %f \n", cimag(alfa));
  printf("iHeight  = %d \n", iHeight);
  printf("PixelWidth  = %f \n", PixelWidth);
  
  printf("distorsion of image  = %f \n", ratio);
  printf("iterMax = %lu \n", iterMax);
  printf("dMaxDistance2Alfa= %f\n", sqrt(dMaxDistance2Alfa2));
  printf("iMaxDistance2Alfa= %d\n", iMaxDistance2Alfa);
 

  // ------------------------------------
  return 0;
}




/* -----------------------------------------  main   -------------------------------------------------------------*/
int main()
{

  
  setup();
  
  


  // here are procedures for creating image file
  
  //FillArray( data ); // no symmetry
  FillArraySymmetric(data);
  SaveArray2PGMFile(data , iterMax+0.000); // save edge array (only boundaries) to pgm file 
  AddBoundaries(data);
  //CheckOrientation( data );
  SaveArray2PGMFile(edge ,iterMax+0.001); // save edge array (only boundaries) to pgm file 
  CopyBoundaries();
  SaveArray2PGMFile(data , iterMax+0.002); // save edge array (only boundaries) to pgm file 
 
  DrawCriticalOrbit(1000000);
  SaveArray2PGMFile(data , iterMax+0.003); // save edge array (only boundaries) to pgm file
  //
  free(data);
  free(edge);
  
  //
  info();

  return 0;
}
