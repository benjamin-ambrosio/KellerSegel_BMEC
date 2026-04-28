#include <iostream>
#include <math.h>
#include <fstream>
#include "Chemo_gmres.h"
#include <ctime>
#include <vector>
using namespace std; 
//This Program computes the solution v of an elliptic problem correspoding to a data driven Keller-Segel equation using the GMRES algorithm

    
        

// Laplacian with Neumann BC (simplified)
double delta(const std::vector<std::vector<double>>& x, int i, int j)
{
    // Mirror indices for Neumann BC
    int ip = std::min(i + 1, nby - 1);
    int im = std::max(i - 1, 0);
    int jp = std::min(j + 1, nbx - 1);
    int jm = std::max(j - 1, 0);

    double sum = x[im][j] + x[ip][j] + x[i][jm] + x[i][jp];
    
    // Subtract coefficient for central point
 
    int neighbors = 4;

    
    return (sum - neighbors * x[i][j]) / (h*h);
}

// Gradient in x-direction with Neumann BC
double gradx(const std::vector<std::vector<double>>& x, int i, int j)
{
    // Mirror indices at boundaries
    int jp = std::min(j + 1, nbx - 1);  // next index
    int jm = std::max(j - 1, 0);        // previous index

    // Approximate gradient using central difference where possible
    return (x[i][jp] - x[i][jm]) / (2.0 * h);
}

// Gradient in y-direction with Neumann BC
double grady(const std::vector<std::vector<double>>& x, int i, int j)
{
    // Mirror indices at boundaries
    int ip = std::min(i + 1, nby - 1);  // next index
    int im = std::max(i - 1, 0);        // previous index

    // Approximate gradient using central difference
    return (x[ip][j] - x[im][j]) / (2.0 * h);
}

   
double operator_AA(const std::vector<std::vector<double>>& x,
                   const std::vector<std::vector<double>>& y, int i, int j)
 { return (x[i][j]*delta(y,i,j)+gradx(x,i,j)*gradx(y,i,j)+grady(x,i,j)*grady(y,i,j));    
   } 

double norm(const std::vector<std::vector<double>>& x)
 { double norm=0;
  for (int i=0;i<nbx;i++)
    for (int j=0;j<nbx;j++)
        norm+=x[i][j]*x[i][j];
   return sqrt(norm);     
   }
 double dot(const std::vector<std::vector<double>>& x,
                   const std::vector<std::vector<double>>& y)
 { double temp=0;
  for (int i=0;i<nbx;i++)
    for (int j=0;j<nbx;j++)
        temp+=x[i][j]*y[i][j];
   return temp;     
   }
 



// basic constructor
syst::syst()
    : t(0),
      u(nbx, std::vector<double>(nby, 0.0)),
      v(nbx, std::vector<double>(nby, 0.0))
{
}



// operator surcharge
syst& syst::operator=(const syst& r)
{
    if (this != &r) // check self-assignment
    {
        t = r.t;
        u = r.u; // std::vector handles copying
        v = r.v;
    }
    return *this;
}





//Initialization operator 
syst& syst::init(ifstream& ff,double bb)
{t=0;
 double aa=0;
  for(int i=0;i<nby;i++)  
    { for(int j=0;j<nbx;j++)
       {ff>>aa;
           cout<<aa<<" ";
       u[i][j]=5*aa;
       v[i][j]=bb;
       }
cout<<endl;}
return *this;
}





//Initialization Operator: file+ function
syst& syst::init(ifstream& ff, double(*pf2)(int k, int l))
{t=0;
 double aa=0;
  for(int i=0;i<nby;i++)  
    { for(int j=0;j<nbx;j++)
       {ff>>aa;
           cout<<aa<<" ";
       u[i][j]=5*aa;// this 5 scaling relates to what was done in recovering file from images
       v[i][j]=(*pf2)(i,j);
       }
cout<<endl;}
return *this;
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


void syst::F(ofstream& output_file )
{ std::cout<<"entering F()";
std::vector<std::vector<std::vector<double>>> V;//contains all the directional vectors vi in the GMRES algorithm
std::vector<std::vector<double>> H(2, std::vector<double>(1, 0.0));// the matrix H in the gmres algorithm, 
std::vector<std::vector<double>> vj(nbx, std::vector<double>(nby, 0.0));// current vector
std::vector<std::vector<double>> B(nbx, std::vector<double>(nby, 0.0));//(1/b)(uf-u-f(u)-d_u\Delta(u))
std::vector<std::vector<double>> r(nbx, std::vector<double>(nby, 0.0));// error
std::vector<std::vector<double>> x(nbx, std::vector<double>(nby, 0.0));// will contain the final v
std::vector<std::vector<double>> cs;// will contain the values c_k,s_k of the rotations matrices
std::vector<double> Qe1(2,0.0);// the rotation matrix of current iteration applied to e1
std::vector<std::vector<double>> R;// will contain QH the product of rotational matrices applied to H
std::vector<double> y;// will contain the coefficients for v=v+y1v1+y2v2+...+ykvk

double nr0;
double aal=0; // in v=v+aal*p
double error=0;
double past_error=0;
double deriv=0;
double deriv2=0;
double normm=0;
double temp=0;
double temp2=0;
double h11=0;
double hkj=0;
double h21=0;
double c1=0;
double s1=0;
double ck=0;
double sk=0;
int count=0;
double vmin;

// we must have B+operator_AA close to 0
// we compute the r which is the error
for(int i=0;i<nby;i++)
    {
      for(int j=0;j<nbx;j++)
        {B[i][j]=(1.0/b)*(ufinal[i][j]-u[i][j]-f(u[i][j])-du*delta(u,i,j)); 
         r[i][j]=operator_AA(u,v,i,j)+B[i][j];
        }
     }

// First step compute v1, h11, h21 and y 

  nr0=norm(r);
error=nr0;
cout<<"Error r0="<<error<<endl;
// optional but recommended
cs.reserve(max_iter);
V.reserve(max_iter);
//H.reserve(max_iter + 1);   // rows



// create V[0]=v1 do not confuse with the v of keller segel
V.emplace_back(nby, std::vector<double>(nbx, 0.0)); // similar yo push_back
for (int i = 0; i < nby; i++)
    for (int j = 0; j < nbx; j++)
        V[0][i][j] = r[i][j] / nr0; // Filling the first Krylov vector
for (int i = 0; i < nby; i++)
    for (int j = 0; j < nbx; j++)
       { h11+=operator_AA(u,V[0],i,j)*V[0][i][j] ;
           cout<<"<Av,v>="<<h11<<endl;}
for (int i = 0; i < nby; i++)
    for (int j = 0; j < nbx; j++)
       { vj[i][j]=operator_AA(u,V[0],i,j)-h11*V[0][i][j] ;
           cout<<"vtilde(i,j)"<<vj[i][j]<<endl;
         }// vtilde_1
h21=norm(vj);

H[0][0]=h11;
H[1][0]=h21; //filling the first column of H

s1=h21/sqrt(h11*h11+h21*h21);
c1=-h11/sqrt(h11*h11+h21*h21);
 Qe1[0]=c1;
Qe1[1]=s1;
cs.emplace_back(std::vector<double>(2, 0.0)); // similar to push_back, here we create v2
cs[0][0]=c1;
cs[0][1]=s1;
V.emplace_back(nby, std::vector<double>(nbx, 0.0)); // similar to push_back, here we create v2
for (int i = 0; i < nby; ++i)
    for (int j = 0; j < nbx; ++j)
        V[1][i][j] = vj[i][j] / h21; // Filling v2 the second Krylov vector

error= nr0*s1;
cout<<"Error r1="<<error<<endl;
count=1;

while((error>0.1)&&(count<max_iter))
{
 //ensure enough rows
H.emplace_back(std::vector<double>(count, 0.0)); // similar to push_back, 

cout<<" We enter the while"<<endl;
for (int l = 0; l < count+1; l++)
{hkj=0;
  for (int i = 0; i < nby; i++)
    for (int j = 0; j < nbx; j++)
        hkj+=operator_AA(u,V[count],i,j)*V[l][i][j] ;
 H[l].push_back(hkj);// 
cout<<"hkj="<<hkj<<endl;
 }          
       
cout<<" First for ok (V)"<<endl;

for (int i = 0; i < nby; i++)
    for (int j = 0; j < nbx; j++)
       { vj[i][j]=operator_AA(u,V[count],i,j);
         for (int l = 0; l < count+1; l++) 
                  vj[i][j]=vj[i][j]-H[l][count]*V[l][i][j];
      
    
         }// vtilde_1
cout<<" Second for ok (Vtilde)"<<endl;
h21=norm(vj);
H[count+1].push_back(h21);
cout<<" H ok "<<endl;
for (int i=0;i<(count+2);i++)
 {for (int j=0;j<(count+1);j++)
    cout<<" H[" <<i<< "][" <<j <<"]="<<H[i][j];
   cout<<endl;}
V.emplace_back(nby, std::vector<double>(nbx, 0.0)); // similar to push_back, here we create V[count+1]

for (int i = 0; i < nby; ++i)
    for (int j = 0; j < nbx; ++j)
        V[count+1][i][j] = vj[i][j] / h21; // Filling v[count+1] the Krylov vector

cout<<" Filling V ok "<<endl;
s1=H[count+1][count]/sqrt(H[count+1][count]*H[count+1][count]+H[count][count]*H[count][count]);
c1=-H[count][count]/sqrt(H[count+1][count]*H[count+1][count]+H[count][count]*H[count][count]);
cs.push_back(std::vector<double>(2, 0.0)); // add the latest to cs (rotation matrix)
cs[count][0]=c1;
cs[count][1]=s1;
cout<<" cs ok "<<endl;
temp=s1*Qe1[count-1]+c1*Qe1[count];
Qe1.push_back(temp);
Qe1[count]=c1*Qe1[count-1]-s1*Qe1[count];// defines the product of rotational matrices to e1; this acts only on the two last rows
cout<<" Qe1 ok "<<endl;
//compute error now


error=abs(Qe1[count+1]);
count++;

cout<<"Error="<<error<<",count="<<count<<endl;

} //end while (that grows the Krylov space)

cout<<" Filling R ok "<<endl;
//REPRENDRE ICI: appliquer les matrices de rotation//REPRENDRE ICI: appliquer les matrices de rotation
R.assign(count, std::vector<double>(count, 0.0));
R[0][0]=sqrt(H[0][0]*H[0][0]+H[0][0]*H[0][0]);
for (int k=1; k<count; k++)// loop on the columns of R
{temp=cs[0][0]*H[0][k]-cs[0][1]*H[1][k];
 R[0][k]=temp;
 cout<<" Filling R[0]["<<k<<"] ok "<<endl;
 for (int i=1;i<k+1;i++)
     {temp=cs[i-1][1]*H[i-1][k]+cs[i-1][0]*H[i][k];
      temp=cs[i][0]*temp-cs[i][1]*H[i+1][k];
      R[i][k]=temp;}
} //end loop k


for (int i=0;i<count;i++)
 {for (int j=0;j<count;j++)
    cout<<" R[" <<i<< "][" <<j <<"]="<<R[i][j];
   cout<<endl;}

y.assign(count,0.0);
//We compute the y, that minimizes the error ||B-Vy|| 
y[count-1]=nr0*Qe1[count-1]/R[count-1][count-1];
for (int i=count-2;i>-1;i--)
  { temp=0;
    for (int j=i+1;j<count;j++)
       temp+=y[j]*R[i][j];
    y[i]=(nr0*Qe1[i]-temp)/R[i][i];
    }

// we compute the new v=x=v+Vy
for (int i = 0; i < nby; i++)
    {for (int j = 0; j < nbx; j++)
       { x[i][j]=v[i][j];
         for (int k=0;k<count;k++)
           x[i][j]+=y[k]*V[k][i][j];
         
         cout<<x[i][j]<<" ";}
      output_file<<endl;
      cout<<endl;}

//Compute the min value of x
vmin=x[0][0];

for (int i = 0; i < nby; i++)
    {for (int j = 0; j < nbx; j++)
       { if (x[i][j]<vmin)
           vmin=x[i][j];  
            }
         }

//Translate x to make it positive
if (vmin<0)
{for (int i = 0; i < nby; i++)
    {for (int j = 0; j < nbx; j++)
       { 
           x[i][j]=x[i][j]-vmin;  
            }
         }
}

// we  write in f the new v=x=v+Vy
for (int i = 0; i < nby; i++)
    {for (int j = 0; j < nbx; j++)
       { output_file<<x[i][j]<<" ";
         cout<<x[i][j]<<" ";}
      output_file<<endl;
      cout<<endl;}

} //end F()

int main()
{

std::cout << "main entered" << std::endl;
   

srand(time(nullptr));


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

ufinal =new double*[nby];
for( int i=0;i<nby;i++)
  {ufinal[i] = new double[nbx];
   for(int j=0;j<nbx;j++)
   {
     ufinal[i][j]=0.0;
    }
  }

ufinal =new double*[nbx];  
for(int i=0;i<nby;i++)
  {ufinal[i] = new double[nbx];
   for(int j=0;j<nbx;j++)
   {
     ufinal[i][j]=0.0;
    }
  }
  std::cout<<"ok3"; 
   
ifstream ff("matrix_Imaget14.dat",ios::in); //fichier recueillant toutes    
ifstream ff2("matrix_Imaget16.dat",ios::in); //fichier recueillant toutes    
ofstream g("ChemoRK4.dat",ios::out); //fichier recueillant toutes 
//les valeurs en espace a des intervalles de temps fixe
ofstream g2("ChemoRK4v.dat",ios::out); //all the values of v 
ofstream g3("EllipSolution_v14.dat",ios::out); //all the values of v 

ofstream h("Chemoev.dat",ios::out);// fichier receuillant les valeurs 
//pour tout temps en certains points d'espace fixes
ofstream k("Chemox.dat",ios::out);

for( int i=0;i<nby;i++)
  {for(int j=0;j<nbx;j++)
   {double a=0.0;
    ff2>>a;
    ufinal[i][j]=5*a;
      }
}


syst S;
//S.init(ff,ff);
 std::cout<<"ok3";
   

   S.init(ff,tau12quar);
//S.init(tau1,tau1bis);
 //while(S.fin())
   //{ 
  
   S.F(g3);
   //if ( (S.div3()))
   //S.printt(h);
  //if (S.div3())
   //S.printux(k);
   //if (S.div3()) // si t est un multiple de 0.1 on ecrit dans g
   S.printu(g);     
  S.printv(g2);
    //}

  g.close();
 g2.close();
  h.close();
   k.close();
   g3.close();
// destruction du tableau de nombres aleatoires 

for( int k=0;k<nby;k++)
{
for( int l=0;l<nbx;l++)
delete[] alea[k][l];
delete[] alea[k];
}   
delete[] alea;
return 0;
}

