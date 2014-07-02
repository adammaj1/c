
/*
  http://www.tenouk.com/Bufferoverflowc/Bufferoverflow1.html
  The output, when the input is: 12345678 (8 bytes), the program run smoothly.
 

  When the input is: 123456789 (9 bytes), then :
  * the following will be displayed when compiled with Microsoft Visual C++ 6.0.  
  * In Linux the “Segmentation fault” message will be displayed and the program terminates.


  gcc b.c -Wall

  a@zalman:~/c/varia/errors/buffer$ ./a.out
  strcpy() NOT executed....
  Syntax: ./a.out <characters>
  a@zalman:~/c/varia/errors/buffer$ ./a.out 1234
  mybuffer content= 1234
  strcpy() executed...
  a@zalman:~/c/varia/errors/buffer$ ./a.out 12345678
  mybuffer content= 12345678
  strcpy() executed...
  a@zalman:~/c/varia/errors/buffer$ ./a.out 123456789
  mybuffer content= 123456789
  strcpy() executed...
  *** stack smashing detected ***: ./a.out terminated
  Przerwane


*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  // theoretically reserve 5 byte of buffer plus the
  // terminating NULL....should allocate 8 bytes = 2 double words,
  // to overflow, need more than 8 bytes...
  // so, if more than 8 characters input by user,
  // there will be access violation, segmentation fault etc.
  char mybuffer[5];
  // a prompt how to execute the program...
  if (argc < 2)
    {
      printf("strcpy() NOT executed....\n");
      printf("Syntax: %s <characters>\n", argv[0]);
      exit(0);
    }
  // copy the user input to mybuffer, without any bound checking
  // a secure version is srtcpy_s()
  strcpy(mybuffer, argv[1]);
  printf("mybuffer content= %s\n", mybuffer);
  // you may want to try strcpy_s()
  printf("strcpy() executed...\n");
  return 0;
}

