/*

gcc -lm -Wall p.c
./a.out


It measures how many iterates needs slowly escaping point z1 
to escape. ( reach target set )

here :
z1 = zp + distance
distance = pow(2.0,-n)

http://www.exploringbinary.com/a-table-of-negative-powers-of-two/

The program works up to n=26
For n=27 it falls into infinite loop.
I think that then it treats z1 as zp
because of numerical problems.


What do you think about it ?

http://arxiv.org/abs/math/0505036
Parabolic Julia Sets are Polynomial Time Computable by Mark Braverman

distance = 2^(-n) and LastIteration is approx 2^n


n=   1 distance = 5.000000e-01  = 0.5000000000 LI =            3 log2(LI) =   2 time =  0 seconds
n=   2 distance = 2.500000e-01  = 0.2500000000 LI =            5 log2(LI) =   2 time =  0 seconds
n=   3 distance = 1.250000e-01  = 0.1250000000 LI =           10 log2(LI) =   3 time =  0 seconds
n=   4 distance = 6.250000e-02  = 0.0625000000 LI =           19 log2(LI) =   4 time =  0 seconds
n=   5 distance = 3.125000e-02  = 0.0312500000 LI =           35 log2(LI) =   5 time =  0 seconds
n=   6 distance = 1.562500e-02  = 0.0156250000 LI =           68 log2(LI) =   6 time =  0 seconds
n=   7 distance = 7.812500e-03  = 0.0078125000 LI =          133 log2(LI) =   7 time =  0 seconds
n=   8 distance = 3.906250e-03  = 0.0039062500 LI =          261 log2(LI) =   8 time =  0 seconds
n=   9 distance = 1.953125e-03  = 0.0019531250 LI =          518 log2(LI) =   9 time =  0 seconds
n=  10 distance = 9.765625e-04  = 0.0009765625 LI =        1 031 log2(LI) =  10 time =  0 seconds
n=  11 distance = 4.882812e-04  = 0.0004882812 LI =        2 055 log2(LI) =  11 time =  0 seconds
n=  12 distance = 2.441406e-04  = 0.0002441406 LI =        4 104 log2(LI) =  12 time =  0 seconds
n=  13 distance = 1.220703e-04  = 0.0001220703 LI =        8 201 log2(LI) =  13 time =  0 seconds
n=  14 distance = 6.103516e-05  = 0.0000610352 LI =       16 394 log2(LI) =  14 time =  0 seconds
n=  15 distance = 3.051758e-05  = 0.0000305176 LI =       32 778 log2(LI) =  15 time =  0 seconds
n=  16 distance = 1.525879e-05  = 0.0000152588 LI =       65 547 log2(LI) =  16 time =  0 seconds
n=  17 distance = 7.629395e-06  = 0.0000076294 LI =      131 084 log2(LI) =  17 time =  0 seconds
n=  18 distance = 3.814697e-06  = 0.0000038147 LI =      262 156 log2(LI) =  18 time =  0 seconds
n=  19 distance = 1.907349e-06  = 0.0000019073 LI =      524 301 log2(LI) =  19 time =  0 seconds
n=  20 distance = 9.536743e-07  = 0.0000009537 LI =    1 048 590 log2(LI) =  20 time =  0 seconds
n=  21 distance = 4.768372e-07  = 0.0000004768 LI =    2 097 167 log2(LI) =  21 time =  0 seconds
n=  22 distance = 2.384186e-07  = 0.0000002384 LI =    4 194 322 log2(LI) =  22 time =  0 seconds
n=  23 distance = 1.192093e-07  = 0.0000001192 LI =    8 388 675 log2(LI) =  23 time =  0 seconds
n=  24 distance = 5.960464e-08  = 0.0000000596 LI =   16 778 821 log2(LI) =  24 time =  0 seconds
n=  25 distance = 2.980232e-08  = 0.0000000298 LI =   33 604 781 log2(LI) =  25 time =  1 seconds
n=  26 distance = 1.490116e-08  = 0.0000000149 LI =   68 563 774 log2(LI) =  26 time =  1 seconds



n=   1 distance(z1,0.5) = 5.000000e-01  = 0.5000000000000 					LastIteration = 3.000000 
n=  10 distance(z1,0.5) = 9.765625e-04  = 0.0009765625000 					LastIteration = 1031.000000 
n= 100 distance(z1,zp)  = 7.888609e-31  = 0.00000000000000000000000000000078886091 

DBL_EPSILON           = 2.220446e-16  = 0.00000000000000022204460492503130808473
C gives the name DBL_EPSILON to the smallest positive number ε such that 1 + ε ≠ 1 to machine precision. 
Since the significant has 52 bits, it’s clear that DBL_EPSILON = 2-52 ≈ 2.2 × 10-16. 
That is why we say a floating point number has between 15 and 16 significant (decimal) figures.

http://www.johndcook.com/blog/2009/04/06/anatomy-of-a-floating-point-number/



Calculated  Machine epsilon: 1.19209E-07
http://en.wikipedia.org/wiki/Machine_epsilon

*/

 
#include <stdio.h>
#include <math.h> // pow
#include <float.h> // DBL_EPSILON 
#include <time.h>


double GiveLastIteration(double Zx0, double Zy0, double Cx, double Cy, int ER2)
 {
  double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
  double Zx = Zx0;
  double Zy = Zy0;
  double distance; // distance(z1,zp)
  double i=0.0;
  Zx2=Zx*Zx;
  Zy2=Zy*Zy;
  while  (Zx2+Zy2<ER2)  /* ER2=ER*ER */
  {
   Zy=2*Zx*Zy + Cy;
   Zx=Zx2-Zy2 + Cx;
   Zx2=Zx*Zx;
   Zy2=Zy*Zy;
   i+=1.0;
  //distance = Zx-0.5; // distance(z1,zp)
  printf("i= %.0f log2(i)= %.0f distance= %.38e = %.38f\n", i,log2(i), Zx, Zx );
  printf("i= %.0f Zx = %.50e\n", i,Zx);
  }
  return i;
 }


int main()

{
  double distance;
  
   int n;
   // parabolic fixed point zp = zpx = zpy*I = 0.5 for c=1/4
   double zpx = 0.5;
   double zpy = 0.0; 
// z1 = z1x + z1y*I = point of exterior of Julia set but near parabolic fixed point zp
   double z1x;
   double z1y = zpy;
   double cx= 0.25;
   double cy= 0.0;
   // Escape Radius ; it defines target set  = { z: abs(z)>ER}
   // all points z in the target set are escaping to infinity
   double ER=2.0;
   double ER2;

   double LastIteration;


   time_t start,end;
   double dif;
  


    ER2= ER*ER;
 
  n=27;

   //for (n =27; n<101; n++)
   {
     time (&start);
     distance = pow(2.0,-n);
     z1x = zpx + distance;
     LastIteration = GiveLastIteration(z1x,z1y, cx,cy,ER2 );
     time (&end);
     dif = difftime (end,start);
     printf("n= %3d distance = %e  = %.10f LI = %10.0f log2(LI) = %3.0f time = %2.0lf seconds\n",n, distance, distance, LastIteration, log2(LastIteration), dif);
   }

 printf("DBL_EPSILON   = %e  = %40.38f ",DBL_EPSILON , DBL_EPSILON );
 

  return 0;
}
