#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// boolean Escape time and DEM/M in one loop  


// gcc t.c -lm -Wall


int main()

{  // limits 
  double ER = 2.0;
  double ER2;
  double absdZMax2= 1e60; 

 // C = Cx + Cy* I = point of parameter c-plane
  double Cx, Cy;
  int i; /* iteration number */
  int iMax = 200;
  unsigned char iInterior = 130;
  unsigned char iExterior = 240; 

  double Zx, Zy; /* Z = Zx + Zy*I  point of dynamical plane */
  double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
  double temp;
  double absZ2;

  // http://en.wikipedia.org/wiki/Complex_quadratic_polynomial
  // first derivative  of fc(zcr) with respect to c =  dZ = dfc(zcr)/dc = 2*Z*dZ  = dZx + dZy*I 
  double dZx = 1.0;
  double dZy = 0.0; 
  double absdZ2; // = abs(dZ)* abs(dZ) = dZx*dZx + dZy*dZy

  double distance;
  unsigned char color = iInterior;

  ER2 = ER * ER;

  Cx=0.25 + 0.0001;
  Cy= 0.0;

  /* initial value of orbit  = critical point zcr = 0 */
  Zx=0.0; 
  Zy=0.0;
  printf(" Zx = %f ; Zy = %f \n", Zx, Zy);
  //
  Zx2 = Zx*Zx;
  Zy2 = Zy*Zy;
  absZ2= Zx*Zx + Zy*Zy;
  absdZ2= dZx*dZx + dZy*dZy;
  // iteration of critical point z= 0 on the dynamical z-plane       
  for (i=0; i<iMax; i++)
    { // check if not escaping : abs(z)>ER
      if (absZ2  > ER2 ) { color= iExterior; break ; } // exterior =  escapes
      if (absdZ2 > absdZMax2) break; // interior when derivative explodes
      // z = fc(z) = z*z + c 
      Zy=2*Zx*Zy + Cy;
      Zx=Zx2-Zy2 +Cx;
      printf("i = %3.0d ; Zx = %f ; Zy = %f \n",i, Zx, Zy);
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

   
   // distance = 2 * |Zn| * log|Zn| / |dZn| 
   // a=sqrt(nz)
   // return 2* a*log(a)/sqrt(nzp) 
   temp = sqrt(absZ2);
   distance = 2.0*temp*log(temp)/sqrt(absdZ2); 
   
   //
   printf(" i = %d ; iMax = %d \n", i, iMax);
   printf(" Zx = %f ; Zy = %f \n", Zx, Zy);
   printf(" color = %u \n", color);
   printf("  true distance = %f and computed distance = %f ; ratio = %f \n", Cx-0.25, distance, (Cx-0.25)/distance);
   printf("  absZ2 = %f ; ER2 = %f \n", absZ2, ER2);
   printf("  absdZ2 = %f ; absdZMax2 = %f \n", absdZ2, absdZMax2);
   printf(" ----------- tests ------------ \n");
   if ( i == iMax ) 
        printf(" ( i == iMax ) = does not escapes or iMax too low  \n ");
        else  printf(" ( i < iMax ) = escapes or derivative explodes  \n ");             
   if (absZ2 > ER2 ) 
        printf(" (absdZ2 > ER2 ) = exterior point which simply escapes ( bailout ) \n ");
        else   printf(" (absdZ2 <= ER2 ) = does not escapes ( bailout ) or iMax too low \n "); 

   if (distance>0.0) 
        printf(" (distance>0.0) = exterior or boundary point \n ");
        else   printf(" (distance<=0.0) = interior point \n ");        
   if (absdZ2 < absZ2) 
       printf(" (absdZ2 < absZ2) = includes escaping through 0\n "); 
       else printf(" (absdZ2 >= absZ2) =  ???\n "); 
   if (absdZ2 > absdZMax2) 
       printf(" (absdZ2 > absdZMax2) =  interior when derivative explodes \n"); // 
       else printf(" (absdZ2 <= absdZMax2) =  ??? \n"); // 

return 0 ; 

}
