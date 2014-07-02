/*

http://stackoverflow.com/questions/575936/how-am-i-incorrectly-using-the-round-function-in-c?rq=1
gcc -lm round.c
./a.out

gcc -lm round.c
./a.out
1.000000
1.000000


gcc 
gcc version 4.8.1 (Ubuntu/Linaro 4.8.1-10ubuntu9) 

*/

#include <stdio.h>
#include <math.h>

int main(void)
{
    float f;
    double d;

    /* Note man page says that roundf() returns a float
       and round() returns a double */
    f = roundf(1.2);
    d = round(1.2);

    printf("%f\n", f);
    printf("%lf\n", d);

    return 0;
}
