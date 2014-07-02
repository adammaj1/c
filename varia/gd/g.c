#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>

gdImagePtr myLoadPng(char *filename)
{
  FILE *in;
  struct stat stat_buf;
  gdImagePtr im;
  in = fopen("myimage.png", "rb");
  if (!in) {
    /* Error */
  }
  if (fstat(fileno(in), &stat_buf) != 0) {
    /* Error */
  }
  /* Read the entire thing into a buffer
    that we allocate */
  char *buffer = malloc(stat_buf.st_size);
  if (!buffer) {
    /* Error */
  }
  if (fread(buffer, 1, stat_buf.st_size, in)
    != stat_buf.st_size)
  {
    /* Error */
  }
  im = gdImageCreateFromPngPtr(
    stat_buf.st_size, buffer);
  /* WE allocated the memory, WE free
    it with our normal free function */
  free(buffer);
  fclose(in);
  return im;
}




