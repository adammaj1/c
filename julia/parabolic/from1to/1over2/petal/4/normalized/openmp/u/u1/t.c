/*

gcc t.c -Wall

 ON THE DERIVATIVE OF THE HAUSDORFF DIMENSION OF THE QUADRATIC JULIA SETS LUDWIK JAKSZTAS
http://www.ams.org/journals/tran/2011-363-10/S0002-9947-2011-05208-6/

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


 Dynamics in one complex variable: introductory lectures by John W. Milnor ( page 7.6 )
http://arxiv.org/abs/math/9201272
normalized : u(zcr) = 0.0

i = 0; z.x = 0.000000; u.x = 0.000000 ; du = NA  
i = 1; z.x = -0.187500; u.x    = 3.120000 ; du = 3.120000  
i = 2; z.x = -0.238998; u.x    = 5.339791 ; du = 2.219791  
i = 3; z.x = -0.269918; u.x    = 7.445030 ; du = 2.105239  
i = 4; z.x = -0.291475; u.x    = 9.498873 ; du = 2.053843  
i = 5; z.x = -0.307719; u.x    = 11.523753 ; du = 2.024880  
i = 6; z.x = -0.320570; u.x    = 13.530317 ; du = 2.006564  
i = 7; z.x = -0.331087; u.x    = 15.524447 ; du = 1.994130  
i = 8; z.x = -0.339912; u.x    = 17.509720 ; du = 1.985273  
i = 9; z.x = -0.347460; u.x    = 19.488463 ; du = 1.978743  
i = 10; z.x = -0.354018; u.x    = 21.462264 ; du = 1.973801  
i = 11; z.x = -0.359786; u.x    = 23.432250 ; du = 1.969987  
i = 12; z.x = -0.364912; u.x    = 25.399247 ; du = 1.966996  
i = 13; z.x = -0.369510; u.x    = 27.363871 ; du = 1.964624  
i = 14; z.x = -0.373664; u.x    = 29.326596 ; du = 1.962725  
i = 15; z.x = -0.377442; u.x    = 31.287790 ; du = 1.961194  
i = 16; z.x = -0.380898; u.x    = 33.247743 ; du = 1.959953  
i = 17; z.x = -0.384076; u.x    = 35.206688 ; du = 1.958945  
i = 18; z.x = -0.387011; u.x    = 37.164812 ; du = 1.958125  
i = 19; z.x = -0.389733; u.x    = 39.122270 ; du = 1.957458  
i = 20; z.x = -0.392266; u.x    = 41.079188 ; du = 1.956918  



*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>



// complex quadaratic polynomial 
double complex f(double complex z, double complex c)
{
  return ( z*z +c);

}


/*

Computing Fatou coordinate with Maxima CAS :
(%i1) z:x+y*%i;
(%o1) %i*y+x
(%i2) a:-1/(2*z^2);
(%o2) −1/(2*(%i*y+x)^2)
(%i3) realpart(a);
(%o3) −x^2−y^2/(2*(y^2+x^2)^2)
(%i4) imagpart(a);
(%o4) (x*y)/(y^2+x^2)^2

*/


double complex GiveU(complex double Z)
{

  double Hx, Hy; 
  double t; // temp 
  // Z = Zx +Zy*I
  double Zx = creal(Z);
  double Zy = cimag(Z);

  // U = Ux +Uy*I
  double Ux; // re(U)  
  double Uy; // im(U) 
 
  // from z to h 
  Hx= Zy;
  Hy= -Zx - 0.5;
  // from h to ux 
  t = Hx*Hx+Hy*Hy;
  Ux = -Hx*Hx - (Hy*Hy)/(2*t*t);
  Uy = (Hx*Hy)/(t*t);
 return (Ux+Uy*I);
  
}

double complex GiveNormalizedU(complex double Z, double ux0)
{
 // U = Ux +Uy*I
  double complex U;
  double Ux; // re(U)  
  double Uy; // im(U) 
 
   U = GiveU(Z);
   Ux = fabs(creal(U)) - fabs(ux0); // re(U)  
   Uy = cimag(U);
   
return (Ux+Uy*I); 


}


int main()
{

 double complex c =  -0.75; 
 double complex z;
 double complex Zcr= 0.0; // critical point
 //double ZAx = -0.5; // alfa fixed point 
 //double dZ;
 //double zx;
 double uxNew, uxOld;
 double ux0;

 double complex z1 = -0.486666666666667  +0.146666666666667*I;  //POints inside Julia set , inside immediate basin of attraction



 int iMax= 20;

int i ; 




z= z1;
ux0 = creal(GiveU(Zcr));
uxOld= creal(GiveNormalizedU(z, ux0));

printf("i = 0; z.x = %f; u.x = %f ; du = NA  \n",  creal(z), uxOld); 

for ( i=1; i<iMax+1; ++i)
{ // 2 iterations 
  z = f(z,c);
  z= f(z,c);
  //
  uxNew = creal(GiveNormalizedU( z, ux0));
  printf("i = %d; z.x = %f; u.x    = %f ; du = %f  \n", i, creal(z), uxNew, (uxNew-uxOld)); 
  uxOld= uxNew; 
} 


return 0;
}
