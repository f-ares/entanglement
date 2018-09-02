/*

 It computes the Rényi entropy of a an interval in 
 the ground state of a quadratic, translational
 invariant Hamiltonian 

*/

/*

 We diagonalise the correlation matrix 
 2<(a_n, a_n^dagger)^t(a_m^dagger, a_m)>-\delta_{nm} 
 using the GSL library

*/

/*

 It reads from a file the entries of a row of this 
 matrix computed with the Mathematica notebook "math_moebius"

*/

/*

 To compile write: 
 gcc -L/usr/local/lib moebius.c -o moebius -lgsl -lgslcblas -lm

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_matrix.h>

FILE *out;

double S, Salpha; /*S is the von Neumann entropy and
                    Salpha is the Rényi entropy of parameter alpha*/

float delta, h,  alpha; /*delta and h are the coupling constants of the theory 
                          and alpha is the Rényi parameter*/

int l; /*Length of the interval*/

void load_correl(gsl_matrix_complex *correl);
void compute_eigenvalues(gsl_matrix_complex *correl, double *eigen);
void compute_entropy(double *eigen);
void write_results();


void main(int argc, char *argv[])
{ 
 if(argc==5)
    {
      sscanf(argv[1], "%d", &l);
      sscanf(argv[2], "%f", &delta);
      sscanf(argv[3], "%f", &h);
      sscanf(argv[4], "%f", &alpha);
    }

 else
    {
      printf("Some parameters are left...\n");
      exit(0);
    }

 gsl_matrix_complex *correl=gsl_matrix_complex_alloc(2*l, 2*l); /*Correlation matrix*/
 double *eigen=malloc(2*l*sizeof(double)); /*Eigenvalues of the correlation matrix*/
  
 load_correl(correl);
 compute_eigenvalues(correl, eigen);
 compute_entropy(eigen); 
 write_results();
      
 free(eigen);
 gsl_matrix_complex_free(correl);
}
   
void load_correl(gsl_matrix_complex *correl) 
{
 int n, m;
 gsl_complex z, w;
 char name_file[256];
 
 
 /*We read the row of the correlation matrix from the file*/

 /*f1 corresponds to the row of the matrix generated by the entries 
     G_00 y G_01 of the symbol*/
 /*f2 corresponds to the row of the matrix generated by the entries 
     G_10 y G_11 of the symbol*/

 double **f1, **f2;

 f1=malloc(2*l*sizeof(int *)); 

 for (n= 0; n < 2*l; n++)     
      f1[n] = malloc(2*sizeof(int));
  
 f2=malloc(2*l*sizeof(int *)); 
  
 for (n= 0; n < 2*l; n++)     
      f2[n] = malloc(2*sizeof(int));
 

 sprintf(name_file, "lrk_correl.dat");
 out=fopen(name_file, "rt");

 for(n=0; n<2*l; n=n+2)
   {
    fscanf(out,"%lf %lf %lf %lf %lf %lf %lf %lf \n", 
      &f1[n][1], &f1[n][2], &f1[n+1][1], &f1[n+1][2], 
      &f2[n][1], &f2[n][2], &f2[n+1][1], &f2[n+1][2]);
   }  
  
 fclose(out);
 
 /*We build the upper triangular part of the correlacion matrix*/ 

 for(n=0; n<2*l; n=n+2)
   {
    for(m=n; m<2*l; m++)
      {
       GSL_SET_COMPLEX(&z, f1[m-n][1], f1[m-n][2]); 
       gsl_matrix_complex_set(correl, n, m, z); 
       GSL_SET_COMPLEX(&z, f2[m-n][1], f2[m-n][2]); 
       gsl_matrix_complex_set(correl, n+1, m, z);   
      }
    }

 /*We build the lower triangular part of the correlation matrix
     taking into account that it is Hermitian*/
    
 for(m=0; m<2*l; m=m+2)
   {
    for(n=m+2; n<2*l; n++)
      {
       gsl_matrix_complex_set(correl, n, m, 
           gsl_complex_conjugate(gsl_matrix_complex_get(correl, m, n)));
       gsl_matrix_complex_set(correl, n, m+1, 
           gsl_complex_conjugate(gsl_matrix_complex_get(correl, m+1, n)));          
      }
   }

 free(f1);
 free(f2);   
}


void compute_eigenvalues(gsl_matrix_complex *correl, double *eigen)
{
 int n;
 gsl_eigen_herm_workspace *w=gsl_eigen_herm_alloc(2*l);
 gsl_vector *eigen2=gsl_vector_alloc(2*l);

 gsl_eigen_herm(correl, eigen2, w);

 for(n=0; n<2*l; n++)
      eigen[n]=gsl_vector_get(eigen2, n);

 gsl_vector_free(eigen2);
 gsl_eigen_herm_free(w);
}

void compute_entropy(double *eigen)
{
 int n;

 S=0;
 Salpha=0;

 for(n=0; n<2*l; n++)
   {
    if(fabs(eigen[n])>=1.0)
      {
       S+=0;
       Salpha+=0;
      }	   
    else
      {
	S+=-(0.5*(1+eigen[n])*log(0.5*(1+eigen[n]))
             +0.5*(1-eigen[n])*log(0.5*(1-eigen[n]))); 

	Salpha+= 1/(1-alpha)*log(pow(0.5*(1+eigen[n]), alpha)
                                 +pow(0.5*(1-eigen[n]), alpha));
      }
    }

  S=0.5*S;
  Salpha= 0.5*Salpha;

  printf("S=%lf\n", S);
}
 
void write_results()
{
  char name_file[256];

  sprintf(name_file, "lrk_entropy_%1.2f_d_%1.2f_h_%1.2f.dat", alpha, delta, h);
  out=fopen(name_file, "at");
  fprintf(out,"%d %.10lf %.10lf \n", l, S, Salpha);
  fclose(out);
}  
