/*

gcc -lm -Wall testf3.c
./a.out


http://arxiv.org/abs/math/0505036
Parabolic Julia Sets are Polynomial Time Computable by Mark Braverman




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



http://randomascii.wordpress.com/2011/08/13/faster-fractals-through-algebra/

*/

 
#include <stdio.h>
#include <math.h> // pow
#include <float.h> // FLT_EPSILON
#include <time.h>
#include <fenv.h> //fegetround()





int main()

{

  
   
   
   float cx= 0.25;
   // Escape Radius ; it defines target set  = { z: abs(z)>ER}
   // all points z in the target set are escaping to infinity
   float ER=2.0;
   float ER2;

   


   time_t start,end;
   float dif;
  


    ER2= ER*ER;
 
   
   
   
   float Zpx = 0.5; 
   float Zx ;// bad value = 0.5002; good value = 0.5004
   float Zx2; /* Zx2=Zx*Zx */
   float i=0.0;
   float d ; // distance between Zpx=1/2 and zx
   float d2; // d2=d*d;
   int n = 13;

   d = pow(2.0,-n);
   Zx = Zpx + d;
   d2=d*d;


    time (&start);
    Zx2=Zpx*Zpx + 2.0*d*Zpx +d2;
    printf ("Using c with float and Escape Radius = %f \n",  ER);
    printf("Round mode is    = %d \n",fegetround());
    printf("i= %3.0f; Zx = %f;  Zx2 = %10.8f ;  d = %f ;  d2 = %.10f\n", i,Zx, Zx2, d , d2);
    if (d2<(FLT_EPSILON/2.0)) 
       {printf("error : relative error to big and d2= %.10f is smaller then FLT_EPSILON/2.0 = %.10f; increase precision ! \a\n", 
                d2, FLT_EPSILON/2.0);    
         return 1;}
  
    while  (Zx2<ER2)  /* ER2=ER*ER */
    {
    
     Zx=Zx2 + cx;
     d = Zx-Zpx;
     d2=d*d; 
     Zx2= 0.25+ d +d2; // zx2 = zx*zx = (zp+d)*(zp+d)= zp2 +2*d*zp +d2 = 2.25+d +d2
     i+=1.0;
     //printf("i= %3.0f; Zx = %f;  Zx2 = %10.8f ;  d = %f ;  d2 = %f \n", i,Zx, Zx2, d,d2);
    }
   
   time (&end);
   dif = difftime (end,start);
   printf("n = %d; distance = %3f; LI = %10.0f log2(LI) = %3.0f time = %2.0lf seconds\n",n, d, i, log2(i), dif);

 
 

  return 0;
}
