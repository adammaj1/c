/*

gcc t.c -Wall



not normalized 
i = 0; z = 0.000000; u = -2.000000  
i = 1; z = -0.025000; u = -2.216066  
i = 2; z = -0.050000; u = -2.469136  
i = 3; z = -0.075000; u = -2.768166  
i = 4; z = -0.100000; u = -3.125000  
i = 5; z = -0.125000; u = -3.555556  
i = 6; z = -0.150000; u = -4.081633  
i = 7; z = -0.175000; u = -4.733728  
i = 8; z = -0.200000; u = -5.555556  
i = 9; z = -0.225000; u = -6.611570  
i = 10; z = -0.250000; u = -8.000000  
i = 11; z = -0.275000; u = -9.876543  
i = 12; z = -0.300000; u = -12.500000  
i = 13; z = -0.325000; u = -16.326531  
i = 14; z = -0.350000; u = -22.222222  
i = 15; z = -0.375000; u = -32.000000  
i = 16; z = -0.400000; u = -50.000000  
i = 17; z = -0.425000; u = -88.888889  
i = 18; z = -0.450000; u = -200.000000  
i = 19; z = -0.475000; u = -800.000000  
i = 20; z = -0.500000; u = -nan  


normalized : u(zcr) = 0.0

i = 0; z = 0.000000; u = 0.000000  
i = 1; z = -0.025000; u = 0.216066  
i = 2; z = -0.050000; u = 0.469136  
i = 3; z = -0.075000; u = 0.768166  
i = 4; z = -0.100000; u = 1.125000  
i = 5; z = -0.125000; u = 1.555556  
i = 6; z = -0.150000; u = 2.081633  
i = 7; z = -0.175000; u = 2.733728  
i = 8; z = -0.200000; u = 3.555556  
i = 9; z = -0.225000; u = 4.611570  
i = 10; z = -0.250000; u = 6.000000  
i = 11; z = -0.275000; u = 7.876543  
i = 12; z = -0.300000; u = 10.500000  
i = 13; z = -0.325000; u = 14.326531  
i = 14; z = -0.350000; u = 20.222222  
i = 15; z = -0.375000; u = 30.000000  
i = 16; z = -0.400000; u = 48.000000  
i = 17; z = -0.425000; u = 86.888889  
i = 18; z = -0.450000; u = 198.000000  
i = 19; z = -0.475000; u = 798.000000  
i = 20; z = -0.500000; u = nan
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double GiveUx(double Zx, double Zy)
{

  double Hx, Hy; 
  double t; // temp 
 // double Ux; // re(fatou coordinate)  
  
  // from z to h 
  Hx= Zy;
  Hy= -Zx - 0.5;
  // from h to ux 
  t = Hx*Hx+Hy*Hy;
  return -Hx*Hx - (Hy*Hy)/(2*t*t);
  
}


double GiveNormalizedUx(double Zx, double Zy, double ux0)
{
 
return fabs(GiveUx(Zx, Zy)) - fabs(ux0);


}


int main()
{

 // c = ( -0.750000 ; 0.000000 ) 
 double Zcrx= 0.0; // critical point
 double ZAx = -0.5; // alfa fixed point 
 double dZ;
 double zx;
 double u;
 double ux0;

 int iMax= 20;

int i ; 




dZ = (Zcrx - ZAx)/(iMax);
ux0 = GiveUx(Zcrx, 0.0);

zx = Zcrx;
for ( i=0; i<iMax+1; ++i)
{  
  zx = Zcrx -  i*dZ;
  u = GiveNormalizedUx( zx, 0.0, ux0);
  printf("i = %d; z = %f; u = %f  \n", i, zx, u); 
   
} 


return 0;
}
