/*

 ix =  0 ; c =  0.25000064588764424551  ; DC =  0.00299935411235575450 ; ratio = 33.34051140812377772918 
 ix =  1 ; c = -0.75000921765452789017  ; DC =  1.00000986354217213567 ; ratio = 0.00299932452839178079 
 ix =  2 ; c = -1.24993551212417045720  ; DC =  0.49992629446964256703 ; ratio = 2.00031459558064224195 
 ix =  3 ; c = -1.36825264215735546529  ; DC =  0.11831713003318500808 ; ratio = 4.22530781746840624832 
 ix =  4 ; c = -1.39401846938445221646  ; DC =  0.02576582722709675117 ; ratio = 4.59201751957555051771 
 ix =  5 ; c = -1.39960549964161574867  ; DC =  0.00558703025716353221 ; ratio = 4.61172144075298967103 
 ix =  6 ; c = -1.40081673950389383227  ; DC =  0.00121123986227808359 ; ratio = 4.61265388562717963726 
 ix =  7 ; c = -1.40107932767567034391  ; DC =  0.00026258817177651165 ; ratio = 4.61269772390573556346 
 ix =  8 ; c = -1.40113625489488162609  ; DC =  0.00005692721921128218 ; ratio = 4.61269978429704084017 
 ix =  9 ; c = -1.40114859630461524413  ; DC =  0.00001234140973361804 ; ratio = 4.61269988113369536036 
 ix = 10 ; c = -1.40115127183305617289  ; DC =  0.00000267552844092876 ; ratio = 4.61269988568461682512 
 ix = 11 ; c = -1.40115185186828653929  ; DC =  0.00000058003523036639 ; ratio = 4.61269988589952526512 
 ix = 12 ; c = -1.40115197761573265830  ; DC =  0.00000012574744611901 ; ratio = 4.61269988590801685059 
 ix = 13 ; c = -1.40115200487687021963  ; DC =  0.00000002726113756133 ; ratio = 4.61269988591366031068 
 ix = 14 ; c = -1.40115201078688783099  ; DC =  0.00000000591001761137 ; ratio = 4.61269988585073266018 
 ix = 15 ; c = -1.40115201206813693824  ; DC =  0.00000000128124910724 ; ratio = 4.61269988634696519788 
 ix = 16 ; c = -1.40115201234590248535  ; DC =  0.00000000027776554711 ; ratio = 4.61269988514696809359 
 ix = 17 ; c = -1.40115201240612004870  ; DC =  0.00000000006021756336 ; ratio = 4.61269987725621798081 
 ix = 18 ; c = -1.40115201241917478113  ; DC =  0.00000000001305473243 ; ratio = 4.61269993101990825358 
 ix = 19 ; c = -1.40115201242200495277  ; DC =  0.00000000000283017164 ; ratio = 4.61269988813091280514 
 ix = 20 ; c = -1.40115201242261851362  ; DC =  0.00000000000061356085 ; ratio = 4.61269919245109646412 
 ix = 21 ; c = -1.40115201242275152918  ; DC =  0.00000000000013301556 ; ratio = 4.61269982035322924016 
 ix = 22 ; c = -1.40115201242278036593  ; DC =  0.00000000000002883674 ; ratio = 4.61271111244792684928 
 ix = 23 ; c = -1.40115201242278661765  ; DC =  0.00000000000000625173 ; ratio = 4.61260448822448059381 
 ix = 24 ; c = -1.40115201242278797291  ; DC =  0.00000000000000135525 ; ratio = 4.61296000000000000008 
 ix = 25 ; c = -1.40115201242278826672  ; DC =  0.00000000000000029382 ; ratio = 4.61254612546125461253 
 ix = 26 ; c = -1.40115201242278833026  ; DC =  0.00000000000000006353 ; ratio = 4.62457337883959044384 
 ix = 27 ; c = -1.40115201242278834414  ; DC =  0.00000000000000001388 ; ratio = 4.57812500000000000000 
 ix = 28 ; c = -1.40115201242278834717  ; DC =  0.00000000000000000304 ; ratio = 4.57142857142857142851 
 ix = 29 ; c = -1.40115201242278834782  ; DC =  0.00000000000000000065 ; ratio = 4.66666666666666666652 



// try to have the same number of the pixels = n
// inside each hyperbolic component of Mandelbrot set along real axis
// width of components 

long double GivePixelWidth(unsigned int period, unsigned int n)
{

  long double w ;
  unsigned int k;

 switch ( period )
 {  // A SCALING CONSTANT EQUAL TO UNITY IN 1D QUADRATIC MAPS M. ROMERA, G. PASTOR and F. MONTOYA
   case      0 : w=(CxMax-CxMin)/n;      break;
   case      1 : w=1.000000000000L/n;    break; // exact value
   case      2 : w=0.310700264133L/n;    break; // numerical approximation , maybe wrong 
   case      4 : w=0.070844843095L/n;    break; // w(2*p) = w(p)/4.0L  ; aproximated value of Feigenbaum constants
   case      8 : w=0.015397875272L/n;    break;
   case     16 : w=0.003307721510L/n;    break;
   case     32 : w=0.000708881730L/n;    break;
   case     64 : w=0.000151841994935L/n; break;
   case    128 : w=0.000032520887170L/n; break;
   case    256 : w=0.00000696502297L/n;  break;
   case    512 : w=0.000001491696694L/n; break;
   case   1024 : w=0.000000319475846L/n; break;
   case   2048 : w=0.000000068421948L/n; break;
   case   4096 : w=0.000000015L/n;       break;
   case   8192 : w=0.000000004L/n;       break;
   case  16384 : w=0.000000001L/n;       break;
   default : if (period == 2*jMax+3)  w=(CxMax-CxMin)/10.0L; // period not found or period > jMax
                else { k=period/16384; w = 0.000000001L; while (k>2) { w /=4.0L; k /=2;};  w /=n;} // dive 
 }

 return w;
}




0.0	0.25            	
1.0	-0.75
2.0	-1.25	
3	-1.3680989
4	-1.3940462
5	-1.3996312
6	-1.4008287			
7	-1.4010853
8	-1.4011402
9	-1.401151982029
10	-1.401154502237
11	-1.401155041988563
12	-1.401155157586859
13	-1.401155182344483
14	-1.40115518764681
15	-1.401155188782406
16	-1.401155189025616
17	-1.401155189077704
18	-1.40115518908886
19	-1.401155189091249
20	-1.401155189091761
21	-1.40115518909187
22	-1.401155189091894
23	-1.401155189091899
24	-1.4011551890919
25	-1.4011551890919
26	-1.4011551890919
27	-1.4011551890919
28	-1.4011551890919
29	-1.4011551890919
30	-1.4011551890919

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



// approximated function using roots on the real axis
long double GiveCx(long double x)
{ // http://zunzun.com
  
  long double a = 5.3226927610784935E-01;
  long double b = 6.5410208763684241E-01;
  long double c = -1.4312869957125389E+00;
  long double d = 8.4710834303177074E-01;
  return (c*atanl(expl((x-a)/b)) + d);
}


long double GiveLyapunov(long double c, unsigned int i, )
{
 long double zx=0.0L; // z = zx+zy*i
 long double zy=0.0L; // critocal point
  


}



int main()
{


long double DC, oldDC;
long double ratio;
long double C, oldC;
int ix;

oldC = 0.25L + 0.003L;
oldDC=0.1L;
ratio = 0.1L;
for (ix=0; ix<30; ix++)
 { 
   C = GiveCx((long double) ix);
   DC = oldC - C;
   ratio = oldDC/DC;
   oldC = C; 
   oldDC = DC;
   printf(" ix = %2d ; c = %.20Lf  ; DC =  %.20Lf ; ratio = %.20Lf \n", ix, C, DC, ratio);
}




return 0;
}
 