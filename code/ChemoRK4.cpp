#include <iostream>
#include <math.h>
#include <fstream>
#include "Chemo.h"
using namespace std; 
// Programme simulant un systeme de type Fitzhugh-Nagumo avec un laplacien et
// un parametre c qui depend du temps et de l'espace. 
//Le systeme est simule grace a la methode des differences finies. 
// La solution est representee par la classe syst ayant pour elements private 
//u,v et t. Les resultats numeriques sont ecrits dans des fichiers.



    
        

// laplacien Neumann
double delta(double** x, int i, int j)
 {
 if (i==0 && j==0)
  return ((x[i+1][j]+x[i][j+1]-2*x[i][j])/(h*h));
 if (i==0 && j== (nbx -1))
  return ((x[i+1][j]+x[i][j-1]-2*x[i][j])/(h*h)); 
 if (i==0 && j!= (nbx -1) && j!=0)
	return ((x[i+1][j]+x[i][j-1]+x[i][j+1]-3*x[i][j])/(h*h));
 if (i==(nby-1)&&j==0)
  return ((x[i-1][j]+x[i][j+1]-2*x[i][j])/(h*h));
 if (i==(nby-1)&& j==(nbx-1))
	return ((x[i-1][j]+x[i][j-1]-2*x[i][j])/(h*h));
 if (i==(nby-1) && j!=0 && j!=(nby-1))
  return ((x[i-1][j]+x[i][j-1]+x[i][j+1]-3*x[i][j])/(h*h));
 if (i!=0 && i!=(nby-1) && j==0 )
  return ((x[i-1][j]+x[i+1][j]+x[i][j+1]-3*x[i][j])/(h*h));
 if (i!=0 && i!=(nby-1) && j==(nbx-1)) 
 return ((x[i-1][j]+x[i+1][j]+x[i][j-1]-3*x[i][j])/(h*h)); 
 if (i!=0 && i!=(nby-1) && j!=0 && j!=(nbx-1)) 
   return ((x[i-1][j]+x[i+1][j]+x[i][j-1]+x[i][j+1]-4*x[i][j])/(h*h));
   }


double gradx(double** x, int i, int j)
 {
 if (j==0)
  return (x[i][j+1]-x[i][j])/h;
 if (j==(nbx -1))
  return (x[i][j]-x[i][j-1])/h; 
 if (j!= (nbx -1) && j!=0)
return 0.5*(x[i][j]-x[i][j-1])/h + 0.5*(x[i][j+1]-x[i][j])/h ; 
   }


double grady(double** x, int i, int j)
 {
 if (i==0)
  return (x[i+1][j]-x[i][j])/h;
 if (i==(nby -1))
  return (x[i][j]-x[i-1][j])/h; 
 if (i!= (nby -1) && i!=0)
return 0.5*(x[i][j]-x[i-1][j])/h + 0.5*(x[i+1][j]-x[i][j])/h ; 
   }




// Constructeur initialise tout a 0
syst::syst()
{  t=0;
 u= new double*[nby];
 v= new double*[nby];
  for(int i=0;i<nby;i++)
    {u[i]= new double[nbx];
     v[i]= new double[nbx];
      for(int j=0;j<nbx;j++)
        {u[i][j]=0;
        v[i][j]=0;



        }
      }
}


// Constructeur
syst::syst(double a,double** x, double** y)
{ t=a;
  u= new double*[nby];
  v= new double*[nby]; 
  for(int i=0;i<nby;i++)
    {u[i]= new double[nbx];
     v[i]= new double[nbx];


      for(int j=0;j<nbx;j++)
        {u[i][j]=x[i][j];
         v[i][j]=y[i][j];
        }
     }
}

//Constructeur
syst::syst(const syst& r)
{t=r.t;
 u= new double*[nby];
  v= new double*[nby];  
 for(int i=0;i<nby;i++)
    {u[i]= new double[nbx];
     v[i]= new double[nbx];

      for(int j=0;j<nbx;j++)
        {u[i][j]=r.u[i][j];
         v[i][j]=r.v[i][j];
         }
      }
}

//surcharge de l'operateur =
syst& syst::operator=(const syst& r)
{t=r.t;
  for(int i=0;i<nby;i++)
  for(int j=0;j<nbx;j++)
   {u[i][j]=r.u[i][j];
    v[i][j]=r.v[i][j];
   }
  return *this;
}

// operateur d'affectation
syst& syst::affect(double a,double** x, double** y)
{t=a;
  for(int i=0;i<nby;i++)
   for(int j=0;j<nbx;j++)
    {u[i][j]=x[i][j];
     v[i][j]=y[i][j];
     }
  return *this;
}

//Second operateur d'affectation: affecte a a t et aux tableaux les valeurs 
//constantes b,c
syst& syst::affect(double aa,double bb, double cc)
{t=aa;
  for(int i=0;i<nby;i++)  
     for(int j=0;j<nbx;j++)
       {u[i][j]=bb;
       v[i][j]=cc;
       }
return *this;
}


//operateur d'initialisation: affecte aux tableaux les valeurs
// constantes b,c
syst& syst::init(double aa,double bb)
{t=0;
  for(int i=0;i<nby;i++)  
     for(int j=0;j<nbx;j++)
       {u[i][j]=aa;
       v[i][j]=bb;
       }
return *this;
}

// second operateur d'initialisation  affecte aux tableaux les valeurs des 
//fonctions *pf1 et *pf2
syst& syst::init(double(*pf1)(int k, int l), double(*pf2)(int k, int l))
{t=0;
  for(int i=0;i<nby;i++)
  for(int j=0;j<nbx;j++)
   {u[i][j]=(*pf1)(i,j);
    v[i][j]=(*pf2)(i,j);
    } 
  return *this;
}

double syst::Moyenneu()
{
 double m=0;
  for(int i=0;i<nby;i++)
    for(int j=0;j<nbx;j++)
       m=m+u[i][j];
  m=m/(nby*nbx);
  return m;
  }
  
double syst::u1kx()
{
 double m=0;
  for(int i=0;i<nby;i++)
    for(int j=0;j<nbx;j++)
       m=m+u[i][j]*cos((double)j*M_PI);
  m=4*m/(nby*nbx);
  return m;
  }
  
double syst::u1ky()
{
 double m=0;
  for(int i=0;i<nby;i++)
    for(int j=0;j<nbx;j++)
       m=m+u[i][j]*cos((double)i*M_PI);
  m=4*m/(nby*nbx);
  return m;
  }
  
double syst::u1kx1ky()
{
 double m=0;
  for(int i=0;i<nby;i++)
    for(int j=0;j<nbx;j++)
       m=m+u[i][j]*cos((double)i*M_PI)*cos((double)j*M_PI);
  m=4*m/(nby*nbx);
  return m;
  }
  
double syst::u2kx()
{
 double m=0;
  for(int i=0;i<nby;i++)
    for(int j=0;j<nbx;j++)
       m=m+u[i][j]*cos(2*(double)j*M_PI);
  m=4*m/(nby*nbx);
  return m;
  }

double syst::u2ky()
{
 double m=0;
  for(int i=0;i<nby;i++)
    for(int j=0;j<nbx;j++)
       m=m+u[i][j]*cos(2*(double)i*M_PI);
  m=4*m/(nby*nbx);
  return m;
  }


double syst::u2kx1ky()
{
 double m=0;
  for(int i=0;i<nby;i++)
    for(int j=0;j<nbx;j++)
       m=m+u[i][j]*cos((double)i*M_PI)*cos(2*(double)j*M_PI);
  m=4*m/(nby*nbx);
  return m;
  }

double syst::u1kx2ky()
{
 double m=0;
  for(int i=0;i<nby;i++)
    for(int j=0;j<nbx;j++)
       m=m+u[i][j]*cos(2*(double)i*M_PI)*cos((double)j*M_PI);
  m=4*m/(nby*nbx);
  return m;
  }
  
double syst::u3ky()
{
 double m=0;
  for(int i=0;i<nby;i++)
    for(int j=0;j<nbx;j++)
       m=m+u[i][j]*cos(3*(double)i*M_PI);
  m=4*m/(nby*nbx);
  return m;
  }

double syst::u3kx()
{
 double m=0;
  for(int i=0;i<nby;i++)
    for(int j=0;j<nbx;j++)
       m=m+u[i][j]*cos(3*(double)j*M_PI);
  m=4*m/(nby*nbx);
  return m;
  }


double syst::Moyennev()
{
 double m=0;
  for(int i=0;i<nby;i++)
    for(int j=0;j<nbx;j++)
       m=m+v[i][j];
  m=m/(nby*nbx);
  return m;
  }

//Ecrit les valeurs de u dans un fichier 
void syst::printu(ofstream& f )
{for(int i=0;i<nby;i++)
   {for(int j=0;j<nbx;j++)
    f<<u[i][j]<<' '; 
   f<< endl; 
   }
f<<endl<<endl;
}

// ecrit v
void syst::printv(ofstream& f )
{for(int i=0;i<nby;i++)
   {for(int j=0;j<nbx;j++)
    f<<v[i][j]<<' '; 
    f<< endl; 
   }
}

// ecrit u et v
void syst::printuv(ofstream& f )
{for(int i=0;i<nby;i++)
   {for(int j=0;j<nbx;j++)
    f<<u[i][j]<<' '; 
   f<< endl; 
   }
 f<<endl<<endl;
for(int i=0;i<nby;i++)
   {for(int j=0;j<nbx;j++)
    f<<v[i][j]<<' '; 
   f<< endl; 
   }
 f<<endl<<endl;
}

// ecrit u(49,j)
void syst::printux(ofstream& f )
{
for (int j=0;j<nbx;j++)
 f<<u[49][j]<<endl;  
  f<<endl;
f<<endl;    
}

void syst::printux2(ofstream& f )
{
for (int j=0;j<nbx;j++)
 f<<u[49][j]<<' ';  
  f<<endl<<endl;
f<<endl;    
}


void syst::printt(ofstream& g )
{
g<<t<<" ";
 for (int i=50;i<80;i++)
   for (int j=50;j<80;j++)
      {g<<u[i][j]<<' '<<v[i][j]<<' '<<delta(u,i,j)<<' '<<delta(v,i,j)<<' '<<gradx(u,i,j)*gradx(v,i,j)+grady(u,i,j)*grady(v,i,j)<<' '<<f(u[i][j]) +du*delta(u,i,j) -b*u[i][j]*delta(v,i,j)   -b*gradx(u,i,j)*gradx(v,i,j)-b*grady(u,i,j)*grady(v,i,j)<<' '<<c*u[i][j]-e*v[i][j]+dv*delta(v,i,j)<<' ';}
      
//g<<Moyenneu()<<' '<<u1kx()<<' '<<u1ky()<<' '<<u2kx()<<' '<<u1kx1ky()<<' '<<u2ky()<<' '<<u3kx()<<' ';
//g<<u2kx1ky()<<' '<<u1kx2ky()<<' '<<u3ky()<<' ';
//f<<u[i_10][j_10]<<' '<<v[i_10][j_10]<<' '<<u[i_4][j_4]<<' '<<v[i_4][j_4];
//f<<' '<<u[i_11][j_11]<<' '<<v[i_11][j_11]<<' ';
 g<<endl;
//if (div2())    //un espace tout les 7, sert   l'affichage avec gnuplot si on 
//veut par exemple u(t). Mais si on veut une animation pour quelques points 
//il faut l'enlever, si on veut affficher tous les temps en une seule fois,
// il faut tout enlever
//f<<endl<<endl;  
}

//destructeur
syst::~syst()
{for(int i=0;i<nby;i++)
  delete[] u[i];
 delete[]u;
for(int i=0;i<nby;i++)
  delete[] v[i];
 delete[]v;
} 
//fonction intervenant dans le calcul de la solution u a l'instant d'apres 

void syst::F()
{
double** x=new double*[nby];
double** y=new double*[nby];

double** k11=new double*[nby];
double** k12=new double*[nby];

double** k21=new double*[nby];
double** k22=new double*[nby];

double** k31=new double*[nby];
double** k32=new double*[nby];

double** k41=new double*[nby];
double** k42=new double*[nby];


for(int i=0;i<nby;i++)
{ 
x[i]=new double[nbx];
y[i]=new double[nbx];
k11[i]=new double[nbx];
k12[i]=new double[nbx];
k21[i]=new double[nbx];
k22[i]=new double[nbx];
k31[i]=new double[nbx];
k32[i]=new double[nbx];
k41[i]=new double[nbx];
k42[i]=new double[nbx];
for(int j=0;j<nbx;j++)
  { k11[i][j]=f(u[i][j]) +du*delta(u,i,j) -b*u[i][j]*delta(v,i,j)   -b*gradx(u,i,j)*gradx(v,i,j)-b*grady(u,i,j)*grady(v,i,j);
    k12[i][j]= c*u[i][j]-e*v[i][j]+dv*delta(v,i,j);
     }
}


for(int i=0;i<nby;i++)
{ 
for(int j=0;j<nbx;j++)
{
k21[i][j]=f(u[i][j]+(dt/2)*k11[i][j])+du*delta(u,i,j)+(dt/2)*du*delta(k11,i,j)
-b*(u[i][j]+(dt/2)*delta(k11,i,j))*(delta(v,i,j)+(dt/2)*delta(k12,i,j)) 
  -b*(gradx(u,i,j)+(dt/2)*gradx(k11,i,j))*(gradx(v,i,j)+(dt/2)*gradx(k12,i,j))
-b*(grady(u,i,j)+(dt/2)*grady(k11,i,j))*(grady(v,i,j)+(dt/2)*grady(k12,i,j));
    
k22[i][j]= c*(u[i][j]+(dt/2)*k11[i][j])-e*(v[i][j]+(dt/2)*k12[i][j])+dv*delta(v,i,j)+dv*(dt/2)*delta(k12,i,j);

   }
}


for(int i=0;i<nby;i++)
{ 
for(int j=0;j<nbx;j++)
{ 
  k31[i][j]=f(u[i][j]+(dt/2)*k21[i][j])+du*delta(u,i,j)+(dt/2)*du*delta(k21,i,j)
-b*(u[i][j]+(dt/2)*delta(k21,i,j))*(delta(v,i,j)+(dt/2)*delta(k22,i,j)) 
  -b*(gradx(u,i,j)+(dt/2)*gradx(k21,i,j))*(gradx(v,i,j)+(dt/2)*gradx(k22,i,j))
-b*(grady(u,i,j)+(dt/2)*grady(k21,i,j))*(grady(v,i,j)+(dt/2)*grady(k22,i,j));
    
k32[i][j]= c*(u[i][j]+(dt/2)*k21[i][j])-e*(v[i][j]+(dt/2)*k22[i][j])+dv*delta(v,i,j)+dv*(dt/2)*delta(k22,i,j);

    }
}


for(int i=0;i<nby;i++)
{ 
for(int j=0;j<nbx;j++)
{ 


k41[i][j]=f(u[i][j]+dt*k31[i][j])+du*delta(u,i,j)+dt*du*delta(k31,i,j)
-b*(u[i][j]+dt*delta(k31,i,j))*(delta(v,i,j)+dt*delta(k32,i,j)) 
  -b*(gradx(u,i,j)+dt*gradx(k31,i,j))*(gradx(v,i,j)+dt*gradx(k32,i,j))
-b*(grady(u,i,j)+dt*grady(k31,i,j))*(grady(v,i,j)+dt*grady(k32,i,j));
    
k42[i][j]= c*(u[i][j]+dt*k31[i][j])-e*(v[i][j]+dt*k32[i][j])+dv*delta(v,i,j)+dt*dv*delta(k32,i,j);


    x[i][j]=u[i][j]+(dt/6)*(k11[i][j]+2*k21[i][j]+2*k31[i][j]+k41[i][j]);
   y[i][j]=v[i][j]+(dt/6)*(k12[i][j]+2*k22[i][j]+2*k32[i][j]+k42[i][j]);
   
}
}



  
cout<<u[nby/2][nbx/2]<<" "<<u[(nby-1)/2][nbx/2]<<" "<<v[1][0]<<" "<<v[(nby-1)/2][nbx/2]<<" "<<delta(u,(nby-1)/2,nbx/2);
cout<<" "<<" "<<t<<endl<<endl;
 affect(t+dt,x,y); 

for (int i=0;i<nby;i++)
{ 
 delete[]k11[i];
   delete[]k12[i];
 delete[]k21[i];
   delete[]k22[i];
 delete[]k31[i];
   delete[]k32[i];
 delete[]k41[i];
   delete[]k42[i]; 
delete[]x[i];
   delete[]y[i];

   }
 delete[]x;
 delete[]y;

 delete[]k11;
 delete[]k12;
delete[]k21;
 delete[]k22;
delete[]k31;
 delete[]k32;
delete[]k41;
 delete[]k42; 
}



  


 int main()
{

   // Generation d'un tableau de nombres "aleatoires" entre 0 et 1

alea =new double**[nby];
for( int k=0;k<nby;k++)
  {alea[k] = new double*[nbx];
   for(int m=0;m<nbx;m++)
   {
   alea[k][m]=new double[2];
       alea[k][m][0]=((double)random()/ (double)(RAND_MAX));
       alea[k][m][1]=((double)random()/ (double)(RAND_MAX));
  }
  }
   
ofstream g("ChemoRK4.dat",ios::out); //fichier recueillant toutes 
//les valeurs en espace a des intervalles de temps fixe
ofstream g2("ChemoRK4v.dat",ios::out); //all the values of v 


ofstream h("Chemoev.dat",ios::out);// fichier receuillant les valeurs 
//pour tout temps en certains points d'espace fixes
ofstream k("Chemox.dat",ios::out);


syst S;
   S.init(tau12bis,tau);
//S.init(tau1,tau1bis);
 while(S.fin())
   { S.F();
   if ( (S.div3()))
   S.printt(h);
  if (S.div3())
   S.printux(k);
   //if (S.div3()) // si t est un multiple de 0.1 on ecrit dans g
   S.printu(g);     
  S.printv(g2);
}

  g.close();
 g2.close();
  h.close();
   k.close();
// destruction du tableau de nombres aleatoires 

for( int k=0;k<nby;k++)
{
for( int l=0;l<nbx;l++)
delete[] alea[k][l];
delete[] alea[k];
}   
delete alea;

}

