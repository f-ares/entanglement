/*

 It computes the logarithm of the determinant of a 
 principal submatrix of the Toeplitz matrix
 obtained in the mathematica notebook "math_matrix_g"

*/

/*To compile this code with GCC:

 gcc -m64 disjoint_mkl.c 
-I/opt/intel/compilers_and_libraries_xxxx.x.xxx/linux/mkl/include  
-L/opt/intel/compilers_and_libraries_xxxx.x.xxx/linux/mkl/lib/intel64_lin 
-L/opt/intel/compilers_and_libraries_xxxx.x.xxx/linux/compiler/lib/intel64 
-lmkl_rt -ldl -lpthread -lm -fopenmp -o disjoint_mkl

*/

/*

In each session before execute the program we must enter in the command line:

export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_xxx.x.xxx
/linux/compiler/lib/intel64:$LD_LIBRARY_PATH 

and

export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2018.0.128/
linux/mkl/lib/intel64:$LD_LIBRARY_PATH  

*/ 

/*

xxxx.x.xxx must be replaced with your version of Intel MKL

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

/*We use the zheev routine from Intel MKL
library in order to diagonalise the correlation
matrix*/

#include <mkl_lapacke.h>
#include <mkl_service.h>

/* Complex datatype */

struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;


FILE *out;

double log_det; 

float nu;

int l1, l2, d, d_max, delta_d;

void load_submatrix(dcomplex* matrix);

extern void zheev(char* jobz, char* uplo, int* n, dcomplex* matrix, int* lda,
                  double* eigen, dcomplex* work, int* lwork, double* rwork, 
                  int* info, int l);

void compute_eigenvalues(dcomplex *submatrix, double *eigen);

void compute_log_det(double autov2[]);

void write_results();

void main(int argc, char *argv[])
{
 int threads;

 if(argc==8)
    {
      sscanf(argv[1],"%d", &l1);
      sscanf(argv[2],"%d", &l2);
      sscanf(argv[3],"%d", &d);
      sscanf(argv[4],"%d", &d_max);
      sscanf(argv[5],"%d", &delta_d);
      sscanf(argv[6],"%f", &nu);
      sscanf(argv[7],"%d", &threads);
    }

  else
    {
      printf("Some parameters are left...\n");
      exit(0);
    }
  
 /*Migration between cores: disabled*/
   setenv("OMP_PROC_BIND", "TRUE", 1); 
 
 /*Dynamic adjustment disabled*/
   mkl_set_dynamic(0); 
 
 /*Number of threads MKL should use*/
   mkl_set_num_threads(threads); 
 
   printf("Number of threads employed: %d\n", mkl_get_max_threads());
  
   while(d<d_max+1)
    {
     printf("d=%d nu=%.2f\n", d, nu);
        
     double *eigen=malloc((l1+l2)*sizeof(double));
     dcomplex* submatrix=(dcomplex*)malloc((l1+l2)*(l1+l2)*sizeof(dcomplex));
    
     load_submatrix(submatrix);
     compute_eigenvalues(submatrix, eigen);   
     compute_log_det(eigen); 
     write_results();
  
     free(eigen);
     free(submatrix);
     
     d+=delta_d;
    }

}
   
void load_submatrix(dcomplex *submatrix) 
 {
   int n, m;
   FILE *out;
   char name_file[256];  
 
   /*f corresponds to a row of the Toeplitz matrix generate in mathematica
     and saved in a textfile*/ 

   dcomplex* f;
   f=(dcomplex*)malloc((l1+d+l2)*sizeof(dcomplex)); 
  
   sprintf(name_file, "matrix_g_%1.2f.dat", nu);

   out=fopen(name_file, "rt");

   for(n=0; n<l1+d+l2; n++)
    {
     fscanf(out,"%*d %lf %lf \n", &f[n].re, &f[n].im); 
    }  
   fclose(out);

   /* * in %*d is works as the assignement-supression character*/
   /*In order to construct the submatrix we must follow the directives 
     of the Intel MKL library*/
  
   /*Block: (1,l1)-(1,l1)*/

   /*We construct the upper triangular part of the block*/ 
   for(n=0; n<l1; n++)
     {
      for(m=n; m<l1; m++)
        {
         submatrix[n*(l1+l2)+m].re=f[m-n].re;
         submatrix[n*(l1+l2)+m].im=f[m-n].im;   
        }
      }

  /*We construct the lower triangular part of the block*/

   for(m=0; m<l1-1; m++)
     {
      for(n=0; n<m+1; n++)
        {
         submatrix[(m+1)*(l1+l2)+n].re=f[m+1-n].re;
         submatrix[(m+1)*(l1+l2)+n].im=-f[m+1-n].im;
        }
     }

  /*Block: (l1+d+1, l1+d+l2)-(1, l1)*/
   for(n=0; n<l2; n++)
    {
     for(m=0; m<l1; m++)
      {
        submatrix[(l1+n)*(l1+l2)+m].re=f[m-(n+l1+d)].re;
        submatrix[(l1+n)*(l1+l2)+m].im=f[m-(n+l1+d)].im;   
       }
     }

  /*Block: (1, l1)-(l1+d+1, l1+d+l2)*/
   for(n=0; n<l1; n++)
    {
     for(m=0; m<l2; m++)
      {
        submatrix[n*(l1+l2)+(m+l1)].re=f[l1+d+m-n].re;
        submatrix[n*(l1+l2)+(m+l1)].im=-f[l1+d+m-n].im;   
       }
     }

  /*Block: (l1+d+1, l1+d+l2)-(l1+d+1, l1+d+l2)*/
  
  /*We construct the upper triangular part of the block*/ 
   for(n=0; n<l2; n++)
     {
     for(m=n; m<l2; m++)
      {
        submatrix[(l1+n)*(l1+l2)+m+l1].re=f[m-n].re;
        submatrix[(l1+n)*(l1+l2)+m+l1].im=f[m-n].im;   
       }
     }

  /*We construct the lower triangular part of the block*/

   for(m=0; m<l2-1; m++)
     {
      for(n=0; n<m+1; n++)
        {
         submatrix[(m+1+l1)*(l1+l2)+n+l1].re=f[m+1-n].re;
         submatrix[(m+1+l1)*(l1+l2)+n+l1].im=-f[m+1-n].im;
        }
     }
 free(f);  
}

void compute_eigenvalues(dcomplex *submatrix, double *eigen)
{
  
  int lwork, l, info;
  dcomplex* work;
  dcomplex wkopt;
  /* rwork dimension should be at least max(1,3*n-2) */
  double *rwork=malloc((3*(l1+l2)-2)*sizeof(double));

  /* Query and allocate the optimal workspace */
        
  lwork = -1;
  l=l1+l2;
  zheev("N", "Lower", &l, submatrix, &l, eigen, &wkopt, &lwork, rwork, &info, l);
  lwork = (int)wkopt.re;
  work = (dcomplex*)malloc(lwork*sizeof(dcomplex) );
        
  /* Solve eigenproblem */
  
  zheev("N", "Lower", &l, submatrix, &l, eigen, work, &lwork, rwork, &info, l);
        
  /* Check for convergence */
  
  if(info > 0) 
    {
     printf("The algorithm failed to compute eigenvalues.\n");
     exit(1);
    }
 
 
  free((void*)work);
  free(rwork);

}

void compute_log_det(double *eigen)
{
  int n;
  
  log_det=0;

  for(n=0; n<l1+l2; n++)
     log_det+=logl(fabsl(eigen[n]));
          
  printf("log_det=%lf\n", log_det);

}
 
void write_results()
{
  char name_file[256];

  sprintf(name_file, "sublog_%.2f_disjoint.dat", nu);
  out=fopen(name_file, "at");
  fprintf(out,"%d %d %d %lf\n", l1, l2, d, log_det);
  fclose(out);  
}


  
