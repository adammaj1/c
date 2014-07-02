/*
c console program

draws dynamic external ray using Newton method 



gcc n.c -lm -Wall
./a.out


*/











#include <stdio.h>
#include <math.h> //  needs -lm also 


// change types from Qt library to standard c 
typedef  unsigned long long int qulonglong; // typedef existing_type new_type_name ;
typedef  unsigned char QColor; // from class type (24 biut :   Rgb or  Hsv ) to 8-bit 
typedef unsigned int uint;

/* 
   http://mndynamics.com/indexm.html 


   qmnplane.cpp by Wolf Jung (C) 2007-2014.
   Defines classes QmnPlane, QmnDraw.

   These classes are part of Mandel 5.10, which is free software; you can
   redistribute and / or modify them under the terms of the GNU General
   Public License as published by the Free Software Foundation; either
   version 3, or (at your option) any later version. In short: there is
   no warranty of any kind; you must redistribute the source code as well.
*/



int rayNewton(//int signtype, 
              uint n, 
              double a, 
              double b,
              //double &x, 
              //double &y, 
              double z[],
              double rlog, 
              double ilog)

{  uint k, l; 
   double fx, fy, px, py, u, v = 0.0, 
   d = 1.0 + z[0]*z[0] + z[1]*z[1], 
   t0, 
   t1;

   for (k = 1; k <= 60; k++)
   {  // if (signtype > 0) { a = x; b = y; } "positive is parameter plane, negative is dynamic plane."
      fx = cos(ilog); 
      fy = sin(ilog);

      t0 = exp(rlog)*fx - 0.5*exp(-rlog)*(a*fx + b*fy);
      t1 = exp(rlog)*fy + 0.5*exp(-rlog)*(a*fy - b*fx);
      fx = z[0]; // x; 
      fy = z[1]; // y; 
      px = 1.0; 
      py = 0.0;

      for (l = 1; l <= n; l++)
      {  u = 2.0*(fx*px - fy*py); 
         py = 2.0*(fx*py + fy*px);
         px = u; 
         // if (signtype > 0) px++;
         u = fx*fx; 
         v = fy*fy; 
         fy = 2.0*fx*fy + b; 
         fx = u - v + a;
         u += v; 
         v = px*px + py*py; 
         if (u + v > 1.0e100) return 1;
      }
      fx -= t0; 
      fy -= t1; 
      if (v < 1.0e-50) return 2;
      u = (fx*px + fy*py)/v; 
      v = (fx*py - fy*px)/v;
      px = u*u + v*v; 
      if (px > 9.0*d) return -1;
      z[0] -= u; 
      z[1] += v; 
      d = px; 
      if (px < 1.0e-28 && k >= 5) break;
   }
   return 0;
} //rayNewton




/* The class mndAngle represents an angle (in turns) by a fraction of
   unsigned long long integers , which are <= 2^64 - 1 . They are normalized
   to the denominator  2^K * (2^P - 1)  when possible,  i.e.,  K + P <= 64 .
   Otherwise  P = 0 . Many applications only need the static functions.
   normalize()   returns the period.
   twice()       doubles the angle modulo 1.
   conjugate()   computes the conjugate angle,  returns period of the char.pt.
   wake()        computes the angles of a limb.
   radians()     gives the angle in radians.
   lambda()      computes  e^{core entropy h}  for dyadic angles.
//   realSpider()  checks if the root is real and computes the center  c .
   truncatedTuning()   truncates  (u1, u2)*w  to the next  n-periodic angle.
*/



//Time ~ nmax^2 , therefore limited  nmax .
// int QmnPlane::newtonRay
// c=a+b*i
int newtonRay(
     //int signtype, 
     qulonglong N, // 
     qulonglong D,
     double a, 
     double b, 
     int q, 
     QColor color, 
     int mode) //5, white, 1

{  uint n; 
   int j, i, k, i0, k0; 
   mndAngle t; 
   
   t.setAngle(N, D, j);
   double logR = 14.0, x, y, u;
   u = exp(0.5*logR); 
   y = t.radians();
   x = u*cos(y); 
   y = u*sin(y);
   if (pointToPos(x, y, i0, k0) > 1) i0 = -10;
   stop(); 
   QPainter *p = new QPainter(buffer); 
   p->setPen(color);

   for (n = 1; n <= (nmax > 5000u ? 5000u : nmax + 2); n++)
   {  t.twice(); // next ray 
      for (j = 1; j <= q; j++)
      {  rayNewton( n, a, b, x, y, exp(-j*0.69315/q)*logR, t.radians()); 
         { n = (n <= 20 ? 65020u : 65010u); break; }
         if (pointToPos(x, y, i, k) > 1) i = -10;
         if (i0 > -10)
         {  if (i > -10) p->drawLine(i0, k0, i, k);
            else { n = 65020u; break; }
         }
         i0 = i; 
         k0 = k;
      }
   }
   //if rayNewton fails after > 20 iterations, endpoint may be accepted
   delete p; 
   update(); 
   if (n >= 65020u) return 2;
   if (mode == 2) { a = x; b = y; }
   if (n >= 65010u) return 1; else return 0;
} //newtonRay








int main () {
 
printf("Hello, world!\n");
 
return 0;
 
}
