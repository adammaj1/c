/*

gcc s.c -Wall
./a.out
*/

#include <stdio.h>
#include <string.h>

int main()
{
    char str1[ 32 ] = "Jam jest ";
    char str2[] = "C.";
    char *comment = "#                                                                                                 ";
    
    printf( "%s\n", strncat( str1, str2, 20 ) );
    printf( "%s\n", strncat( comment, str2, 2 ) );
    return 0;
}
