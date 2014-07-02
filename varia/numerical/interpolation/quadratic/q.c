
// gcc q.c -Wall
// ./a.out

# include <stdio.h>



// input values : 3 complex points in 2 arrays 
double xin[3]={3,1,-2}; // double xin[3]={x0,x1,x2};
double yin[3]={10,0,15}; // double yin[3]={y0,y1,y2};


double DeterminatOfMatrix33(double w[3][3])
{
return ( w[0][0]*w[1][1]*w[2][2] + w[1][0]*w[2][1]*w[0][2] + w[2][0]*w[0][1]*w[1][2] - w[2][0]*w[1][1]*w[0][2] - w[0][0]*w[2][1]*w[1][2] - w[1][0]*w[0][1]*w[2][2] );
  
}


double GiveDet_n(int n, double ws[3][3], double wy[3])
{
  int i;
  double wi[3][3]; // use local copy, do not change ws !
  
 
 // copy values from ws to wi
 for (i=0; i<3; ++i)
   {
   wi[i][0]= ws[i][0];
   wi[i][1]= ws[i][1];      
   wi[i][2]= ws[i][2];		
  }
  
  // copy wy column
 for (i=0; i<3; ++i)
   wi[i][n]=wy[i];
   
 return DeterminatOfMatrix33(wi);
}


// main matrix of system of equations 
double GiveMatrixOfSystem(double wx[3], double ws[3][3])
{
 int i;
 
 //printf(" ws = {");
 for (i=0; i<3; ++i)
  {
   ws[i][0]= wx[i]*wx[i]; //printf("{ %f ,",ws[i][0]);
   ws[i][1]= wx[i];       //printf("%f ,",ws[i][1]);
   ws[i][2]= 1;		  //printf("%f }, \n",ws[i][2]);
  }
  
  

return DeterminatOfMatrix33(ws);
}

// find coefficients ( out) of quadratic functions 
// which goes thru 3 point ( xin and yin )
// result is in out array
// uses matrix method of solving system of 3 equations ( linear)   
int FindFunction(double xin[3], double yin[3], double out[3] )
{

 double ws[3][3];
 double dets,deta, detb, detc;

 // compute dets and ws 
 dets = GiveMatrixOfSystem(xin,ws);
 // compute det? 
 deta = GiveDet_n(0,ws,yin);
 detb = GiveDet_n(1,ws,yin);
 detc = GiveDet_n(2,ws,yin);
 
 // compute coeff
 out[0] = deta/dets;
 out[1] = detb/dets;
 out[2] = detc/dets;


return 0;
}

/* =================== main ============================================================*/

int main()
{
 double out[3];

 FindFunction(xin, yin, out );
 // print input
 printf("z0 = ( %f ; %f ); \n",xin[0], yin[0]);
 printf("z1 = ( %f ; %f ); \n",xin[1], yin[1]);
 printf("z2 = ( %f ; %f ); \n",xin[2], yin[2]);
 // print output : coefficients of quadratic function thru 3 points
 printf("a = %f ; b = %f ; c = %f ;\n",out[0], out[1], out[2]);
 
 return 0;
}


