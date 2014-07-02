

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// approiximated function using roots on the real axis
long double GiveCx(long double x)
{ // http://zunzun.com
  
  long double a = 0.53226927610784935L;
  long double b = 0.65410208763684241L;
  long double c = -1.4312869957125389L;
  long double d = 0.84710834303177074L;
  return (c*atanl(expl(x-a)/b) + d);
}

double BurkardtCollectionBased_sech_cdf_Offset_model(double x_in)
{
double temp;
temp = 0.0;
// coefficients
double a = 5.3226927610784935E-01;
double b = 6.5410208763684241E-01;
double c = -1.4312869957125389E+00;
double Offset = 8.4710834303177074E-01;
temp = c * atan(exp((x_in-a)/b));
temp += Offset;
return temp;
}

int main()
{




int ix;

for (ix=0; ix<5; ix++)
printf(" ix = %d ; c = %.20f  ;  %.20Lf \n", ix, BurkardtCollectionBased_sech_cdf_Offset_model( (long double)ix), GiveCx((long double )ix));




return 0;
}
 
