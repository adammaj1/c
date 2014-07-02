#include <stdio.h>
 
int main(){
  int iX,iY;
  const int iXmax = 300; 
  const int iYmax = 300;
  /* color  is coded from 0 to 255 */
  /* it is 8 bit color RGB file */
  const int MaxColorComponentValue=255; 
  FILE * fp;
  char *filename="m.pgm";
  char *comment="# this is my new text pgm file ";  /* comment should start with # */
  char *comment2=" end of comment ";  /* comment should start with # */
  static unsigned char color;
 
 
 
  /*create new file,give it a name and open it in text mode  */
  fp= fopen(filename,"w"); /*  text mode */
  /*write ASCII header to the file*/
  fprintf(fp,"P2\n%s %s\n%d %d\n%d\n",comment,comment2, iXmax,iYmax,MaxColorComponentValue);
  /*write image data bytes to the file*/
  for(iY=0;iY<iYmax;++iY){
    for(iX=0;iX<iXmax;++iX){         
      color=150;   /* compute  pixel color (8 bit = 1 byte) */
      fprintf(fp," %d ", color);   /*write color to the file*/
    }
    fprintf(fp," \n ");
  }
  fclose(fp);
  printf("OK\n");
  return 0;
}
