#ifndef FITZ_H
#define FITZ_H
#include <iostream>
#include <math.h>
#include <fstream>
#include <time.h>
#include <stdlib.h>
// This will be seen by the compiler only once 

using namespace std;
const double T=50 ; // Temps maximal unite= min

const double h=1;//pas d'espace unite = 0.1 mm

const double du=1.0, dv=10.0;// du=1 good value
const double b=10.0, a=7, e=2.0, c=3, gam=0.25; // al=10 good value, gam=0.29 gives BMEC patterns 

const int nbx=100; // taille de la longueur du tableau (solution en espace)
const int nby=100;
 //constantes
double dt=0.01; // pas de temps unite = 1 min
double dtmin=dt;

double*** alea;
double** alea2;













// classe représentant la solution a  un instant t donne
class syst {
public:
  syst();
  syst(double, double **, double **);
  syst(const syst&);
  syst& operator=(const syst&);
  syst& affect(double, double **, double **);
  syst& affect(double, double, double);
  syst& init(double, double);
  syst& init(double (*)(int,int), double (*)(int,int)); 
  void F();
  void EulerImplicite();
  void RemplitM(double**,double*);//sert à remplir la matrice et le vecteur intervenant dans l'algorithme de Newton pour Euler Implicite
void RemplitNLM(double**,double*);
void RemplitV(double*,double*,double);
void RemplitNLV(double*,double*,double);
void printu(ofstream& );
void printv(ofstream& );
  void printuv(ofstream& );
void printux(ofstream& );
void printux2(ofstream& );
  void printt(ofstream& );
  int fin(){return t<T;};
  int grand(){return t>180;};
 int cent(){return ((int)t==20);};
 int cent2(){return ((int)t==1090);};
  int div(){return (floor(10*t)-floor(10*(t-dt)));};
  int div2(){return ((int)(100*t)% 700)==0;};
int div3(){return ((int)(1000*t)% 100 )==0;};
int div4(){return (floor(100*t)-floor(100*(t-dt)));};
double Moyenneu();
double u1kx();// first coefficient in the eigenbasis of Laplacian in the square
double u1ky();// idem....
double u1kx1ky();// idem....
double u2kx();
double u2ky();
double u3kx();
double u3ky();
double u2kx1ky();
double u1kx2ky();
double Moyennev();
double Moyenneeq1();
double Moyenneeq2();
double Normeeq1();
double Normeeq2();
  ~syst();
private:
 double t;  
 double** u;
 double** v;       
};











// fonction servant initialiser
double initi(int i, int j)
{if ((i==49)&&(j==49))
return 0;
else if ((i==49)&&(j==50))
return 1;
else if ((i==50)&&(j==50))
return 2;
else if ((i==50)&&(j==49))
return 3;
else return 0.5;
}




//fonction cubique 
double f(double x)
{return (double) a*x*(x-gam)*(1.0-x);} //GOOD FOR BMEC PATTERNS with gam=0.29

// fonction cubique 
//double f(double x)
//{return (double) a*x*(x-0.08)*(1.0-x);}

// fonction abs value
double vabs(double x)
{ if (x>0)
   return x;
  else 
   return -x;}

double Ber(double p)
{double u=((double)random()/ (double)(RAND_MAX));
if (u<p)
return 1;
else 
return 0;
}

double Ber(int i,int j, double p)
{
if (alea[i][j][0]<p)
return 1;
else 
return 0;
}

double Bin(double p, int n)
{double s,b;
for (int i=1;i<n;i++)
{b=Ber(p);
s+=b;
}
return s;
}

double Bin(int i,int j, double p, int n)
{
return Bin(p,n);
}


double Normcr()
{double s;
srand((unsigned)time(0) );
for (int i=1;i<24;i++)
{
s+=((double)rand()/ (double)(RAND_MAX));
}
s=(s-12)/(12*sqrt(2));
return s;
}


double ALEAT(int i,int j)
{
return (Ber(i,j,0.5)-0.5)*0.01*(0.1+alea[i][j][0]*0.9);
}

double epsilon1(int i,int j)
{if (Ber(i,j,0.5)==1) 
return 0.1;
else 
return 0.15;
}

double epsilon2(int i,int j)
{if ((i>50|i==50)&&(((double)(i-50)>(-sqrt(3)*(j-50)))|(j==50)))
    return 0.1;
  else if ((i<50)&&(((double)(i-50)<(j-50)*sqrt(3))|(j==50)))
    return 0.15 ;
  else 
    return 0.2;
}

double epsilon3(int i,int j)
{return 0.02+ALEAT(i,j);
}
double epsilon(int i,int j)
{
if ((i<nby/2-1)|(i>nby/2)|(j<nbx/2-1)|(j>(nbx/2)))
    
      return 0.015;
else
	return 0.025;
   }



double Norme2(double* x, int N)
{double S=0;
 for(int i=0; i<N;i++)
  S=S+x[i]*x[i];
return sqrt(S);
 }



//fonction donnant la valeur du parametre c en fonction du temps et de l'espace
double c0(double t, int i, int j)
{
return 0;
  }

//fonction donnant la valeur du parametre c en fonction du temps et de l'espace
double c1(double t, int i, int j)
{return -1.5;
  }






//fonction donnant la valeur du parametre c en fonction du temps et de l'espace
double c3(int i, int j)
{
if (j==0)
 return 1;
else return 0;
 }




double c7(double t, int i, int j)
{// nombre aleatoires entre 0 et 2
   return 2*(alea[i][j][0]); 
   
}

double tau1(int i, int j)
{
   return 0.2*cos((1.0/120.0)*M_PI*(double)j)*cos((1.0/80.0)*M_PI*(double)i); 
   
}

double tau1bis(int i, int j)
{
   return 0.8*cos((1.0/120.0)*M_PI*(double)j)*cos((1.0/80.0)*M_PI*(double)i); 
   

}


double tau7(int i, int j)
{
if ((i>10)&&(i<20))
    {if ((j>10)&&(j<20))
      return 0.8;
    else if ((j>40)&&(j<50))
      return 0.8;
    else if ((j>70)&&(j<80))
      return 0.8;
      else 
	return 0.2;}
else if ((i>40)&&(i<50))
  {if ((j>10)&&(j<20))
      return 0.8;
    else if ((j>40)&&(j<50))
      return 0.8;
    else if ((j>70)&&(j<80))
      return 0.8;
      else 
	return 0.2;//0.2
   }

else if ((i>70)&&(i<80))
  {if ((j>10)&&(j<20))
      return 0.8;
    else if ((j>40)&&(j<50))
      return 0.8;
    else if ((j>70)&&(j<80))
      return 0.8;
      else 
	return 0.2;}
else  
    return 0.2;
}




double tau7bis(int i, int j)
{
if ((i>10)&&(i<20))
    {if ((j>10)&&(j<20))
      return 0.8;
    else if ((j>30)&&(j<40))
      return 0.8;
    else if ((j>50)&&(j<60))
      return 0.8;
     else if ((j>70)&&(j<80))
      return 0.8;
      else 
	return 0.0;}
else if ((i>30)&&(i<40))
   {if ((j>10)&&(j<20))
      return 0.8;
    else if ((j>30)&&(j<40))
      return 0.8;
    else if ((j>50)&&(j<60))
      return 0.8;
     else if ((j>70)&&(j<80))
      return 0.8;
      else 
	return 0.0;}

else if ((i>50)&&(i<60))
   {if ((j>10)&&(j<20))
      return 0.8;
    else if ((j>30)&&(j<40))
      return 0.8;
    else if ((j>50)&&(j<60))
      return 0.8;
     else if ((j>70)&&(j<80))
      return 0.8;
      else 
	return 0.0;}

else if ((i>70)&&(i<80))
  {if ((j>10)&&(j<20))
      return 0.8;
    else if ((j>30)&&(j<40))
      return 0.8;
    else if ((j>50)&&(j<60))
      return 0.8;
     else if ((j>70)&&(j<80))
      return 0.8;
      else 
	return 0.0;}
else  
    return 0.0;
}





double tau8(int i, int j)
{
if (i<20)
    {if (j<20) 
      return 0.2;
    else if (j<40)
      return 0.5;
    else if (j<60)
      return 0.2;
      else 
	return 0.5;}
else if (i<40)
  {if (j<20) 
      return 0.5;
    else if (j<40)
      return 0.2;

    else if (j<60)
      return 0.5;
      else 
	return 0.2;}

else if (i<60)
{if (j<20) 
     return 0.2;
    else if (j<40)
      return 0.5;
    else if (j<60)
      return 0.2;
      else 
	return 0.5;}
else if (i<80)
  {if (j<20) 
      return 0.5;
    else if (j<40)
      return 0.2;

    else if (j<60)
      return 0.5;
      else 
	return 0.2;}
else 
{if (j<20) 
     return 0.2;
    else if (j<40)
      return 0.5;
    else if (j<60)
      return 0.2;
      else 
	return 0.5;}
}


double tau9(int i, int j)
//4 ondes spirales
{
if (i<25)
    {if (j<25) 
      return 0.2;
    else if (j<50)
      return -1;
    else if (j<75)
      return -1;
      else 
	return 0;}
else if (i<50)
  {if (j<25) 
      return 1;
    else if (j<50)
      return 0;

    else if (j<75)
      return 0;
      else 
	return 1;}

else if (i<75)
{if (j<25) 
      return 1;
    else if (j<50)
      return 0;
    else if (j<75)
      return 0;
      else 
   return 1;}	
else 
{if (j<25) 
      return 0;
    else if (j<50)
      return -1;
    else if (j<75)
      return -1;
      else 
   return 0;}	
}
double tau9bis(int i, int j)
{
if (i<25)
    {if (j<25) 
      return -1;
    else if (j<50)
      return 0;
    else if (j<75)
      return 0;
      else 
	return -1;}
else if (i<50)
  {if (j<25) 
      return 0;
    else if (j<50)
      return 1;

    else if (j<75)
      return 1;
      else 
	return 0;}

else if (i<75)
{if (j<25) 
      return 0;
    else if (j<50)
      return 1;
    else if (j<75)
      return 1;
      else 
   return 0;}	
else 
{if (j<25) 
      return -1;
    else if (j<50)
      return 0;
    else if (j<75)
      return 0;
      else 
   return -1;}	
}
double tau10(int i, int j)
{

if (i<50)
    
      return 1.18;
else

      return -1.18;
    
   
}

double tau10bis(int i, int j)
{

if (i<50)
    
      return -0.4;
else

      return  0.4;
    
   
}




double tau12(int i, int j)
{
// simulation d'une loi exponentielle

  return (-50)*log(1-alea[i][j][0]);

}

double tau12ter(int i, int j)
{


  return 1-2*alea[i][j][1];

}

double tau12bis(int i, int j)
{


  //return 1.0+ 0.005*alea[i][j][0];
return tau7(i,j) + 0.2*alea[i][j][0]; //Good for BMEC

}


double tau12p5(int i, int j)
{
// simulation d'une loi normale centree reduite

  return sqrt((-2.0)*log(alea[i][j][0]))*cos(2*M_PI*alea[i][j][1]);

}

double tau12p6(int i, int j)
{
// simulation d'une loi normale NON centree reduite

  return 0.5+sqrt((-2.0)*log(alea[i][j][0]))*cos(2*M_PI*alea[i][j][1]);

}
double tau13(int i, int j)
{
// on obtient une spirale 
  if ((i>50|i==50)&&(((double)(i-50)>(-sqrt(3)*(j-50)))|(j==50)))
    return 1.7;
  else if ((i<50)&&(((double)(i-50)<(j-50)*sqrt(3))|(j==50)))
    return -2.02 ;
  else 
    return -1.1;
}

double tau13bis(int i, int j)
{ 
   if ((i>50)&&(((double)(i-50)>(-sqrt(3)*(j-50)))|(j==50)))
    return 0.13;
   else if ((i<50|i==50)&&(((double)(i-50)<(j-50)*sqrt(3))|(j==50)))
    return 2.22;
  else  
    return  -1.8;
}

double tau14(int i, int j)
{// train d'ondes
 
if (j<(nbx/2))
    return -2.5;
  else
    return 2.5;

}

double tau14bis(int i, int j)
{

   if (j<(nbx/2))
    return 2.5;
   else 
    return -2.5;
}

double tau15(int i, int j)
{
// on obtient une spirale dans le cas oscillant
  if ((i>50|i==50)&&(((double)(i-50)>(-sqrt(3)*(j-50)))|(j==50)))
    return 1.18;
  else if ((i<50)&&(((double)(i-50)<(j-50)*sqrt(3))|(j==50)))
    return 1 ;
  else 
    return -0.82;
}

double tau15bis(int i, int j)
{ // on obtient une spirale dans le cas oscillant
   if ((i>50)&&(((double)(i-50)>(-sqrt(3)*(j-50)))|(j==50)))
    return -0.6;
   else if ((i<50|i==50)&&(((double)(i-50)<(j-50)*sqrt(3))|(j==50)))
    return 0.6;
  else  
    return  -0.46;
}

double tau15ter(int i, int j)
{
// on obtient une spirale dans le cas oscillant
  if ((i>50|i==50)&&(((double)(i-50)>(-sqrt(3)*(j-50)))))
    return 0;
  else if ((i<50)&&(((double)(i-50)<(j-50)*sqrt(3))))
    return (a+b)/100 ;
  else 
    return (a+b)/50;
}


double tau16(int i, int j)//1 spirale
{
// 
  if (i<(nbx/2))
    {if (j<(nbx/2))
        return 1;
     else
         return 0;
        }
 else
    {if (j<(nbx/2))
        return 0;
    
     else
         return -1;
        }
}

double tau16bis(int i, int j)
{ 
  if (i<(nbx/2))
    {if (j<(nbx/2))
        return 0;
     else
         return 1;
        }
 else
    {if (j<(nbx/2))
        return -1;
    
     else
         return 0;
        }
}


//fonction donnant l'avance de c(t) sur le chemin suivant la position d'espace
double tau(int i, int j)
{return 0.5;
}









/*

FONCTIONS SERVANT DANS EULER IMPLICITE
*/





















#endif
