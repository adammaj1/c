/*

gcc -lm -Wall ps2.c
./a.out


http://arxiv.org/abs/math/0505036
Parabolic Julia Sets are Polynomial Time Computable by Mark Braverman


It measures how many iterates needs slowly escaping point z1 
to escape. ( reach target set )

here :
z1 = zp + distance
distance = pow(2.0,-n)

http://www.exploringbinary.com/a-table-of-negative-powers-of-two/

A floating-point number is composed of four elements:

    a sign: either negative or non-negative
    a base (or radix): which expresses the amount of quantities that can be represented with a single digit (2 for binary, 10 for decimal, 16 for hexadecimal, and so on...)
    a significand (or mantissa): which is a series of digits of the abovementioned base. The number of digits in this series is what is known as precision.
    an exponent (also knon as characteristic, or scale): which represents the offset of the significand, which affects the value in the following way:
    value of floating-point = significand x baseexponent, with its corresponding sign.





distance = 2^(-n) and LastIteration is approx 2^n

Using c with float and Escape Radius = 2.000000 
n=   1 distance = 5.000000e-01  = 0.5000000000 LI =          3 log2(LI) =   2 time =  0 seconds
n=   2 distance = 2.500000e-01  = 0.2500000000 LI =          5 log2(LI) =   2 time =  0 seconds
n=   3 distance = 1.250000e-01  = 0.1250000000 LI =         10 log2(LI) =   3 time =  0 seconds
n=   4 distance = 6.250000e-02  = 0.0625000000 LI =         19 log2(LI) =   4 time =  0 seconds
n=   5 distance = 3.125000e-02  = 0.0312500000 LI =         35 log2(LI) =   5 time =  0 seconds
n=   6 distance = 1.562500e-02  = 0.0156250000 LI =         68 log2(LI) =   6 time =  0 seconds
n=   7 distance = 7.812500e-03  = 0.0078125000 LI =        133 log2(LI) =   7 time =  0 seconds
n=   8 distance = 3.906250e-03  = 0.0039062500 LI =        261 log2(LI) =   8 time =  0 seconds
n=   9 distance = 1.953125e-03  = 0.0019531250 LI =        518 log2(LI) =   9 time =  0 seconds
n=  10 distance = 9.765625e-04  = 0.0009765625 LI =       1032 log2(LI) =  10 time =  0 seconds
n=  11 distance = 4.882812e-04  = 0.0004882812 LI =       2068 log2(LI) =  11 time =  0 seconds
n=  12 distance = 2.441406e-04  = 0.0002441406 LI =       4058 log2(LI) =  12 time =  0 seconds



http://www.cplusplus.com/reference/clibrary/cfloat/
FLT_EPSILON = 1E-5 = Difference between 1 and the least value greater than 1 that is representable.


C gives the name FLT_EPSILON to the smallest positive number ε such that 1 + ε ≠ 1 to machine precision. 
Since the significant has 24 bits, it’s clear that FLT_EPSILON = 1E-5. 
This gives from 6 to 9 significant decimal digits precision (if a decimal string with at most 6 significant decimal is converted to IEEE 754 single precision and then converted back to the same number of significant decimal, then the final string should match the original; and if an IEEE 754 single precision is converted to a decimal string with at least 9 significant decimal and then converted back to single, then the final number must match the original ).

http://www.johndcook.com/blog/2009/04/06/anatomy-of-a-floating-point-number/



Calculated Machine epsilon using float : 1.19209E-07
http://en.wikipedia.org/wiki/Machine_epsilon




"Using higher-precision arithmetic doesn't eliminate numerical instability, it just makes it take longer to become apparent." 
http://stackoverflow.com/questions/8995030/dealing-with-lack-of-floating-point-precision-in-opencl-particle-system

http://forums.nvidia.com/index.php?showtopic=73067

http://www.codeproject.com/Articles/25294/Avoiding-Overflow-Underflow-and-Loss-of-Precision

http://www.johndcook.com/blog/2008/04/16/overflow-and-loss-of-precision/

http://randomascii.wordpress.com/2012/01/23/stupid-float-tricks-2/#more-400
*/

 
#include <stdio.h>
#include <math.h> // pow
#include <float.h> // FLT_EPSILON
#include <time.h>
#include <fenv.h> //fegetround()



float GiveLastIteration(float Zx0, float Zy0, float Cx, float Cy, int ER2)
 {
  float Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
  float Zx = Zx0;
  //float Zy = Zy0;
 // float distance; // distance(z1,zp)
  float i=0.0;
  Zx2=Zx*Zx;
  printf("i= %.0f; Zx = %f;  Zx2 = %f \n", i,Zx, Zx2);
  //Zy2=Zy*Zy;
  while  (Zx2<ER2)  /* ER2=ER*ER */
  {
   //Zy=2*Zx*Zy + Cy;
   Zx=Zx2 + Cx;
   Zx2=Zx*Zx;
   //Zy2=Zy*Zy;
   i+=1.0;
  //distance = Zx-0.5; // distance(z1,zp)
 // printf("i= %.0f log2(i)= %.0f distance= %.38e = %.38f\n", i,log2(i), Zx, Zx );
  printf("i= %.0f; Zx = %f;  Zx2 = %f ;  Cx = %f\n", i,Zx, Zx2, Cx);
  }
  return i;
 }


int main()

{

  
  float distance;
  
   int n;
   // parabolic fixed point zp = zpx = zpy*I = 0.5 for c=1/4
   float zpx = 0.5;
   float zpy = 0.0; 
// z1 = z1x + z1y*I = point of exterior of Julia set but near parabolic fixed point zp
   float z1x;
   float z1y = zpy;
   float cx= 0.25;
   float cy= 0.0;
   // Escape Radius ; it defines target set  = { z: abs(z)>ER}
   // all points z in the target set are escaping to infinity
   float ER=2.0;
   float ER2;

   float LastIteration;


   time_t start,end;
   float dif;
  


    ER2= ER*ER;
 
   
   printf ("Using c with float and Escape Radius = %f \n",  ER);
   printf("FLT_EPSILON   = %e  = %40.38f \n",FLT_EPSILON , FLT_EPSILON );
   printf("Round mode is    = %d \n",fegetround());
   n=12;
   
   //for (n =12; n<101; n++)
   {
     time (&start);
     distance = pow(2.0,-n);
     z1x = 0.50020; //zpx + distance;
     printf("n= %d; Zpx = %f;  z1x = %f ;  z1x2 = %f\n", n,zpx, z1x,  z1x*z1x);
     LastIteration = GiveLastIteration(z1x,z1y, cx,cy,ER2 );
     time (&end);
     dif = difftime (end,start);
     printf("n= %3d distance = %e  = %.10f LI = %10.0f log2(LI) = %3.0f time = %2.0lf seconds\n",n, distance, distance, LastIteration, log2(LastIteration), dif);
   }

 
 

  return 0;
}
