/* 
 * image.c
 *
 * read in a PPM or PGM binary image and display it, full size
 
  Based on code by : 
   Dr. Lori L. Scarlatos
   Stony Brook University
   http://ms.cc.sunysb.edu/~lscarlatos/
   "I do not have a license for image.c; 
   it was created as an example for my students.
   Please feel free to use it.
   Best regards,
   Lori"
* ----------------------------
* it does not opens asci versions of these files 
* examples files : 
* http://people.sc.fsu.edu/~jburkardt/data/data.html

  // gcc v.c -lm -lGLU -lglut -Wall // ? windows 
    gcc v.c -lm -lglut -lGL -lGLU -Wall // ubuntu 
  ./a.out 5.pgm

 */

//#include <Windows.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h> /* fabs*/
//#include <malloc.h>



/* Global variables */

#define MAXLINE 80	/* maximum length of a line of text */


GLint ImageWidth, ImageHeight;	         /* size of the Image in pixels */
GLint WindowWidth, WindowHeight;	/* size of the window in pixels */
GLint ScreenWidth, ScreenHeight;	/* size of the screen in pixels */


GLubyte *Picture;	/* Array of colors (GLubyte)  */
int size; 



// mouse position as a global variables 
//static float mx=0.0f,my=0.0f ;
int iX, iY; // 
int centerX = 200, centerY = 200;
GLint iYmax,iXmax; // mouse coordinate inside image 



// change it manually !!!!!


const double  ZyMin=-1.0;
const double  ZxMin=-2.0;
const double  PixelHeight=0.0020010005002501  ;
const double  PixelWidth=0.0020010005002501  ;

int filetype;
enum {P2, P3, P5, P6};	/* possible file types */

/* gives position of point (iX,iY) in 1D array  ; uses also global variables */
unsigned int f(unsigned int _iX, unsigned int _iY)
{return (_iX + (iYmax-_iY-1)*iXmax );}

/* 
* Read from a PPM or PGM (binary) file 
* Output is a array of GLubyte 
*/
void readPPM (char *filename, GLubyte **pic) {

	FILE *fp;
	char line[MAXLINE];
	int i, rowsize;
        // int size; // moved to global var 
        int j ;
	GLubyte *ptr;

/* Read in file type */
  perror( NULL );
  fp = fopen(filename, "r"); /* in Unix rb = r */
  if (fp==NULL) { printf("Error from fopen : I can't open %s file' ! ", filename); exit(1); }
  else printf("File %s has been opened !\n", filename);
  /* Each file starts with aa two-byte magic number (in ASCII) that explains :
   * - the type of file it is (PBM, PGM, and PPM) 
   * - its encoding (ASCII or binary). 
   * The magic number is a capital P followed by a single digit number. 
  */
  fgets (line, MAXLINE, fp); /* 1st line : Magic Number  */
  switch (line[1])
  {
   case '2':
       filetype = P2; 
       printf("This is PGM text file (P2), but now I do not have procedure for opening it \n");
       break;
   case '3' :
      filetype = P3;
      printf("This is PPM text file (P3), but now I do not have procedure for opening it !\n");
      break; 
   case '5':
       filetype = P5; 
       printf("This is PGM binary file (P5) and I can open it !\n");
       break;
   case '6' :
      filetype = P6;
      printf("This is PPM binary file (P6) and I can open it !\n");
      break;
   default : 
      printf("Error from readPPM : need binary PPM or binary PGM file as input!\n");
      exit(1);
   }
 /* if this is a comment, read next line. Maybe in binary files there is no comment ?*/
/* there maybe more then one line of comment */
 fgets (line, MAXLINE, fp); 
 while  (line[0]=='#') 
  { printf(" comment  = %s \n", line); /* 2nd or more line  starting with # =  comment, print it   */
    fgets (line, MAXLINE, fp); // read next line 
  } 

/* Read in width and height, & allocate space */
/* these 2 numbers should be in one line with space between them */
   /* 3nd line: width  and height */
  sscanf(line, "%d %d", &ImageWidth, &ImageHeight);
  printf ("iWidth = %d\n", ImageWidth);
  printf ("iHeight = %d\n",  ImageHeight);
  iXmax=ImageWidth-1;
  iYmax=ImageHeight-1;

  if (filetype == P5) {
	  size = ImageHeight * ImageWidth; /* greymap: 1 byte per pixel */
	  rowsize = ImageWidth;
  }
  else /* filetype == P6 */ {
	  size = ImageHeight * ImageWidth * 3; /* pixmap: 3 bytes per pixel */
	  rowsize = ImageWidth * 3;
  }
  *pic = (GLubyte *)malloc (size); /* create dynamic array */

/* Read in maximum value (ignore) */
  fgets (line, MAXLINE, fp); /* next  line */
  /*  */
  if (filetype==P5 || filetype==P6){
    /* Read in the pixel array row-by-row: 1st row = top scanline */
    ptr = *pic + (ImageHeight-1) * rowsize;
    for (i = ImageHeight; i > 0; i--) {
          /* For binary File I/O you use fread and fwrite */
	  j = fread((void *)ptr, 1, rowsize, fp); 
	  ptr -= rowsize;
    }
   if (j) printf("File %s has been read !\n", filename);
          else printf(" j Error from readPPM procedure : I can't read %d file !.\n", j);
  }
  else printf("Error from readPPM procedure : I can't read %s file !. It should be P5 or P6 file !\n", filename);
  fclose(fp);
  printf("File %s has been closed !\n", filename);
}






/* Draw the picture on the screen */

void Draw(void) {
        /* black background of GLUT window */
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Clear the background of our window to black
        glClear(GL_COLOR_BUFFER_BIT); //Clear the colour buffer (more buffers later on)
        glLoadIdentity(); // Load the Identity Matrix to reset our drawing locations
        glFlush(); // Flush the OpenGL buffers to the window
        // left lower corner of displayed image 
	glRasterPos2i(-1, -1); // By default, OpenGL assumes a system running from -1 to 1, 
        switch (filetype){
          case P5 : 	/* greymap: use as illumination values */
		glDrawPixels(ImageWidth, ImageHeight, GL_LUMINANCE, GL_UNSIGNED_BYTE, Picture);
                printf("P5 Image has been drawn !\n");
                break;
	  case  P6 :
		glDrawPixels(ImageWidth, ImageHeight, GL_RGB, GL_UNSIGNED_BYTE, Picture);
                printf("P6 Image has been drawn !\n");
                break;
         default : 
                printf("Error from Draw procedure : There is no image to draw !\n");
          }
}




// Detecting Mouse Clicks
// x and y specify the location (in window-relative coordinates) of the mouse when
// the event occurred
void MouseClicks (int button, int state, int x, int y)
{
switch (button)
case GLUT_LEFT_BUTTON:
if (state == GLUT_DOWN)
// do something
case GLUT_RIGHT_BUTTON: ;
// etc., etc.
}


double GiveUx(double Zx, double Zy)
{

  double Hx, Hy; 
  double t; // temp 
 // double Ux; // re(fatou coordinate)  
  
  // from z to h 
  Hx= Zy;
  Hy= -Zx - 0.5;
  // from h to u 
  t = Hx*Hx+Hy*Hy;
  return fabs(-Hx*Hx - (Hy*Hy)/(2*t*t));
  
}


// mouse motion 
//The x and y callback parameters 
// indicate the  mouse  location 
// in window relative coordinates
//
// x, y –> coordinates of the mouse relative to upper left corner of window

// setting the mouse position to be relative to the mouse
// position inside the window
//

void PassiveMouseMotion(int x, int y) {
  
   //double Zx,Zy;
   //double Ux; 
   //int index;
   //GLubyte Gray;

   // put your code here  ???? 
   iX = x;
   /* invert y axis */
   iY = glutGet(GLUT_WINDOW_HEIGHT) - y -1 ; //+ (glutGet(GLUT_WINDOW_HEIGHT) - ImageHeight); ///glutGet(GLUT_WINDOW_HEIGHT) - y; // ????
   
  

   
   /* output : prints to console  but image can be smaller/greater then window */
   if ((filetype==P5 || filetype==P6) && -1<iX && iX< ImageWidth  && iY<ImageHeight) // && iY<ImageHeight
   printf(" coordinates of the mouse relative to upper left corner of window x=%3d  y=%3d  windowHeight = %3d \n", iX, iY, (glutGet(GLUT_WINDOW_HEIGHT)));
}




static void Key(unsigned char key, int x, int y)
{
    switch (key) {
	case 27:  /* esc */
        case 'q':
        case 'Q':{printf("Key preseed and exit \n"); exit(1) ;}
	default: return ;
    }
}



/* Resize the picture  ; OpenGL calls this function whenever the window is resized */
void Reshape(GLint w, GLint h) {
    /* the viewport is the rectangular region of the window where the image is drawn */
   // Viewport : A rectangular region in the screen for display  (in screen coordinate system) 
    glViewport(0, 0, ImageWidth-1, ImageHeight-1); // glViewport( 0.f, 0.f, SCREEN_WIDTH, SCREEN_HEIGHT );
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // Define a world window : A rectangular region in the world that is to be displayed (in world coordinate system)
    // By default, OpenGL assumes a system running from -1 to 1, 
    gluOrtho2D(-1, 1, -1, 1); // An orthographic projection is basically a 3D projection that does not have perspective
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}





/* 
Initialization: create window 
glutInitWindowSize(600, 600);
*/
void MyInit(int argc, char *argv[]) {
    char filename[MAXLINE];
  
  /* Read in the file (allocates space for Picture) */
  if (argc < 2) 
    {
	printf ("Enter the name of a binary PPM or PGM file: ");
	scanf("%s", filename);
	readPPM ((char *)filename, &Picture);
    }
    else { readPPM (argv[1], &Picture); }

    glutInit(&argc, argv);
      
  glutInitWindowPosition(0, 0); // upper left  corner
  glutInitWindowSize(2560, 1440); // full screen of my monitor
  glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
  if (glutCreateWindow("OpenGl Image viewer ") == GL_FALSE) 
      {printf("Error from  MyInit , glutCreateWindow\n"); exit(1);}
    
}


/* ------------------ Main program ---------------------------------------*/
int main(int argc, char **argv)
{
  
    MyInit(argc, argv);
    glutPassiveMotionFunc(PassiveMouseMotion); // mouse move with no key pressed 
    glutReshapeFunc(Reshape); //  move or resize of screen window 
    glutDisplayFunc(Draw); // 
    glutKeyboardFunc(Key);

    glutMainLoop();
    return 0;
}
