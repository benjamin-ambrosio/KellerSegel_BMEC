// Ce programme utilise glut et openGL pour afficher l'évolution en 2D de 
//la solution du probleme resolu dans fitz.cpp

//g++ -framework GLUT -framework OpenGL -framework Cocoa anifitz.cpp
//g++ anifitz.cpp -lGL -lGLU -lglut souslinux mais attention mettre les bons chemins gl.h dans le include


#include </usr/include/GL/gl.h>
#include </usr/include/GL/glu.h>
#include </usr/include/GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>



using namespace std;
/*	Create checkerboard image	*/
#define	L 300
#define	H 300
GLubyte Image[L][H][3];

static GLdouble zoomFactor = 1.0;
static GLint height;
const int nby=100;
const int nbx=100;
ifstream f("ChemoRK4.dat",ios::in);
double t=0.0;



// Fonction qui copie les donnees de finale.dat 
//(que l'on a recuperees dans u**) dans l'image et on grossit  
//en utilisant glpixelzoom
void makeImage(double **u)
{int i, j;
    
   for (i = 0; i < nby; i++) 
      for (j = 0; j < nbx; j++)
         {
          if (u[i][j]<0)// si u<0 noir
  Image[nby-1-i][j][0]= Image[nby-1-i][j][1]= Image[nby-1-i][j][2]= 0;
else if (u[i][j]<1) // bleu
   {Image[nby-1-i][j][0]= Image[nby-1-i][j][1]=0;
    Image[nby-1-i][j][2]=(int) (((u[i][j])/1.0)*255);} 
 else if (u[i][j]<2) //magenta
{Image[nby-1-i][j][0]= (int) (((u[i][j]-1)/1.0)*255);
Image[nby-1-i][j][1]= 0;
 Image[nby-1-i][j][2]= 255;} 
 else if (u[i][j]<3)//rouge
   {Image[nby-1-i][j][0]= 255;
Image[nby-1-i][j][1]= 0;
Image[nby-1-i][j][2]=  (int)(3-u[i][j])*(255/1.0);} 
 else if (u[i][j]<4)//jaune
   {Image[nby-1-i][j][0]=  255;
Image[nby-1-i][j][1]= (int)((u[i][j]-3.0)*(255/1.0));
Image[nby-1-i][j][2]= 0;}
 else
   {Image[nby-1-i][j][0]= Image[nby-i][j][1]= 255;
Image[nby-1-i][j][2]= 0;}
    }
}

void init(void)
{  glClearColor (0.0,0.0,0.0,0.0);
   glShadeModel(GL_FLAT);
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
   

}


void display(void)
{  // On recupere les donees de finale.dat
double ** u;
u= new double*[nby];
for (int i=0;i<nby;i++)
   {u[i]=new double [nbx];
   for (int j=0;j<nbx;j++)
     f>>u[i][j];
   }
makeImage(u);
   glDrawBuffer(GL_BACK);
   glClear(GL_COLOR_BUFFER_BIT);
   glRasterPos2f(0.0,0.0);
   glBitmap(0,0,0,0,(-L)/2,(-H)/2,0);
   glPixelZoom(3,3);//on grossit l'image 3*3
   glDrawPixels(L, H, GL_RGB, GL_UNSIGNED_BYTE, Image);
   cout <<t<<endl;
   glutSwapBuffers();
   glFlush();
t=t+0.1;// On ecrit le temps correspondant a accorder avec fitz.cpp
for (int i=0;i<nby;i++)
delete[]u[i];
delete[]u;
}

void reshape(int w, int h)
{
   glViewport(0, 0, (GLsizei) w, (GLsizei) h);
   height = (GLint) h;
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluOrtho2D(0.0, (GLdouble) w, 0.0, (GLdouble) h);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}


void idle(void) 
{
 glutPostRedisplay();
         }

int main(int argc, char** argv)
{  glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE| GLUT_RGB);
   glutInitWindowSize(L, H);
   glutInitWindowPosition(500,500);
   glutCreateWindow("hello");
   glutDisplayFunc(display);
   init();
   glutIdleFunc(idle);// appelée quand rien ne se passe.
   glutMainLoop();
   return 0; 
}




