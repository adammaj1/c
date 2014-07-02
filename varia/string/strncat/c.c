#include <stdio.h>
#include <string.h>
 
int main()
{
  unsigned int i=100;
  unsigned int length = 30;
  char filename [length];
  
// snprintf(buf, sizeof buf, "%s%s%s", "String1", "String2", "String3");
 // int snprintf(char *str, size_t size, const char *format, ...)
//char str[LEN];
  snprintf(filename, length, "%d%s", i, ".pgm");
  printf(" %s\n", filename);

 FILE * fp;
  fp = fopen(filename,"wb");
  fclose(fp);

  return 0;

}
