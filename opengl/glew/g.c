/*


gcc g.c  -lGLEW -lGL -lGLU -lglut 
./a.out

http://www.linuxforu.com/2013/04/graphics-programming-in-linux/

Graphics Programming in Linux By Rajnish-Singh on April 1, 2013 in Developers, Features, Overview · 3 Comments


http://www.opengl.org/sdk/libs/GLEW/
call glewInit(); once you've obtained a OpenGL context in your program.
http://steps3d.narod.ru/tutorials/glew-tutorial.html
http://www.packtpub.com/article/using-glew

*/

#include <stdio.h>
#include <GL/glew.h> // glew before gl 
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h> 
 


void display() 
{ 
    glClearColor(1,0,0,0); 
    glClear(GL_COLOR_BUFFER_BIT); 
    glFlush(); 
} 


int main(int argc, char**argv) { 
    glutInit(&argc, argv); 
    glutInitWindowPosition(100,100); 
    glutInitWindowSize(500,500); 
    glutCreateWindow("Hello World");

    GLenum err = glewInit();
    if (GLEW_OK != err)
      {
       fprintf(stderr, "Error: %s\n", glewGetErrorString(err));/* Problem: glewInit failed, something is seriously wrong. */
       return 1; 
      }    
    fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
    
    // checking for extensions 
    if ( !GLEW_ARB_shading_language_100 )
    {
        printf ( "GL_ARB_shading_language_100 NOT supported.\n" );

        return 2;
    }
    printf ( "GL_ARB_shading_language_100 is supported.\n" );

    if ( !GLEW_ARB_shader_objects )
    {
        printf ( "GL_ARB_shader_objects NOT supported\n" );

        return 3;
    }
    printf ( "GL_ARB_shader_objects is supported\n" );

    if (!GLEW_ARB_vertex_shader || !GLEW_ARB_fragment_shader)
     {
           printf("No GLSL support\n");
          exit(4);
     }

   printf(" GLSL is supported\n");

 
    glutDisplayFunc(display); 
    glutMainLoop(); 
}
