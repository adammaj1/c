#include <stdio.h>


void TestArray(int length)
{

 double n[length]; /* n is an array of 10 integers */
   int i;
 
   /* initialize elements of array n to 0 */         
   for ( i = 0; i < length; i++ )
   {
      n[ i ] = i *100.0; /* set element at location i to i + 100 */
   }

   n[1]=n[1]; // to remove warning: variable ‘n’ set but not used [-Wunused-but-set-variable

   // http://stackoverflow.com/questions/37538/how-do-i-determine-the-size-of-my-array-in-c
   // length = sizeof(n)/sizeof(n[0])
   printf(" length = %d , size = %ld bytes = %ld kB  \n",length, length*sizeof(double), length*sizeof(double)/1000);
}
   
 


 
int main ()
{
  int i;
for (i=1; i<1000000; i++)
  
  TestArray(1000*i);
   


return 0;
}
