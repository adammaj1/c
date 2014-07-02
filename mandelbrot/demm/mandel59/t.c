#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// boolean Escape time and DEM/M in one loop  


// gcc t.c -lm -Wall
// ./a.out

/*

true distance = 0.003000 and computed distance = 7.743709 ; ratio = 0.000387 
true distance = 0.001500 and computed distance = 7.715447 ; ratio = 0.000194
true distance = 0.000300 and computed distance = 7.692883 ; ratio = 0.000039 

*/



  /* global variables */
   // limits 
  double PixelSize = 0.003;
  double ER = 2.0;
  double ER2;
  double absdZMax2= 1e60; 
 int iMax = 200;
   // point c on 
   double Cxb= -2.0;
   double Cyb = 0.0;
 

  
 

// 
double GiveDistance(double Cx, double Cy, double Cxb, double Cyb)
{  double dx, dy;
   dx= Cx-Cxb;
   dy = Cy-Cyb; 
   return sqrt(dx*dx + dy*dy);
} 


void PrintInfo(double Cx, double Cy, double Cxb, double Cyb)
{

  int i; /* iteration number */
  double Zx, Zy; /* Z = Zx + Zy*I  point of dynamical plane */
  double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
  double temp;
  double absZ2;

  // http://en.wikipedia.org/wiki/Complex_quadratic_polynomial
  // first derivative  of fc(zcr) with respect to c =  dZ = dfc(zcr)/dc = 2*Z*dZ  = dZx + dZy*I 
  double dZx = 0.0;
  double dZy = 0.0; 
  double absdZ2; // = abs(dZ)*abs(dZ) = dZx*dZx + dZy*dZy
  double TrueDistance;
  double EstimDistance;
  

  

  /* initial value of orbit  = critical point zcr = 0 */
  Zx=0.0; 
  Zy=0.0;
  //printf(" Zx = %f ; Zy = %f \n", Zx, Zy);
  //
  Zx2 = Zx*Zx;
  Zy2 = Zy*Zy;
  absZ2= Zx*Zx + Zy*Zy;
  absdZ2= dZx*dZx + dZy*dZy;

   // 
 
  // iteration of critical point z= 0 on the dynamical z-plane       
  for (i=0; i<iMax; i++)
    { // check if not escaping : abs(z)>ER
      if (absZ2  > ER2 )   break ;  // exterior =  escapes
      if (absdZ2 > absdZMax2) break; // interior when derivative explodes
      // z = fc(z) = z*z + c 
      Zy=2*Zx*Zy + Cy;
      Zx=Zx2-Zy2 + Cx;
      printf("i = %d ; Zx = %f ; Zy = %f \n",i, Zx, Zy);
     /* first derivative   dz = 2*z*dz + 1.0 = dZx + dZy*I  */
      temp = 2*(Zx*dZx - Zy*dZy) + 1.0 ;  
      dZy = 2*(Zx*dZy + Zy*dZx); 
      dZx = temp;

      // abs 
      Zx2 = Zx*Zx;
      Zy2 = Zy*Zy;
      absZ2= Zx2 + Zy2; // nz = x*x + y*y;
      absdZ2= dZx*dZx + dZy*dZy; // 
    };

   // x = log(nz); x = x*x*nz / (temp[1]*nzp); //4*square of dist/pixelwidth
   temp = log(absZ2);
   EstimDistance = temp*temp*absZ2/absdZ2; 
   TrueDistance = GiveDistance(Cx, Cy, Cxb, Cyb) ;
   
 //
   printf(" i = %d ; iMax = %d \n", i, iMax);
   printf(" last : Zx = %f ; Zy = %f \n", Zx, Zy);
   printf("  absZ2 = %f ; ER2 = %f \n", absZ2, ER2);
   printf("  absdZ2 = %f ; ratio absdZ2/absdZMax2 = %e \n", absdZ2, absdZ2/absdZMax2);
   printf(" sqrt(absdZ2) = %.16f \n", sqrt(absdZ2));
   printf("  temp = sqrt(absZ2) =  %.16f ; log(temp) = %.16f \n",temp, log(temp));
   printf("  true distance = %f and computed distance = %f ; ratio = %f \n", TrueDistance, EstimDistance, TrueDistance / EstimDistance);
   printf(" ----------- tests ------------ \n");
   if ( i == iMax ) 
        printf(" ( i == iMax ) = does not escapes or iMax too low  \n ");
        else  printf(" ( i < iMax ) = escapes or derivative explodes  \n ");             
   if (absZ2 > ER2 ) 
        printf(" (absdZ2 > ER2 ) = exterior point which simply escapes ( bailout ) \n ");
        else   printf(" (absdZ2 <= ER2 ) = does not escapes ( bailout ) or iMax too low \n "); 

   if (EstimDistance > 0.0) 
        printf(" (EstimDistance > 0.0) = exterior or boundary point \n ");
        else   printf(" (EstimDistance <= 0.0) = interior point \n ");        
   if (absdZ2 < absZ2) 
       printf(" (absdZ2 < absZ2) = includes escaping through 0\n "); 
       else printf(" (absdZ2 >= absZ2) =  ???\n "); 
   if (absdZ2 > absdZMax2) 
       printf(" (absdZ2 > absdZMax2) =  interior when derivative explodes \n"); // 
       else printf(" (absdZ2 <= absdZMax2) =  derivative not explodes  \n"); // 

}


int main()

{  

   

 
  
 
 
  
   double Cx, Cy ;  // parameter c of fc(z) = point on the parameter c-plane
   //int i;   
  ER2 = ER * ER;
  Cx= Cxb - PixelSize*1.0;
  Cy= Cyb;  
  PrintInfo(Cx, Cy, Cxb, Cyb );   

return 0 ; 

}
