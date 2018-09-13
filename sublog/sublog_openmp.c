/*

 It computes the logarithm of the determinant of the Toeplitz matrix
 obtained in the mathematica notebook "math_matrix_g"

*/

/*To compile this code with GCC:

 gcc -m64 sublog_memento_openmp.c 
-I/opt/intel/compilers_and_libraries_xxxx.x.xxx/linux/mkl/include  
-L/opt/intel/compilers_and_libraries_xxxx.x.xxx/linux/mkl/lib/intel64_lin 
-L/opt/intel/compilers_and_libraries_xxxx.x.xxx/linux/compiler/lib/intel64 
-liomp5 -lmkl_rt -ldl -lpthread -lm -fopenmp -o sublog_memento_openmp

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
 
/*We use the zheev routine from Intel MKL
library in order to diagonalise the correlation
matrix*/

#include <mkl_lapacke.h>
#include <mkl_service.h>
#include <omp.h>

/* Intel MKL complex datatype */
struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;


FILE *out;

double log_det; 

float nu;

void diagonalize(int l);

void load_matrix(dcomplex* matrix, int l);

extern void zheev(char* jobz, char* uplo, int* n, dcomplex* a, int* lda,
                  double* eigen, dcomplex* work, int* lwork, double* rwork, 
                  int* info, int l);

void compute_eigenvalues(dcomplex *matrix, double *eigen, int l);

void compute_log_det(double *eigen , int l);

void write_results(int l);

void main(int argc, char *argv[])
{
 int i, delta_l, arg, l, NUM_THREADS;
 
 /*We simultaneously consider NUM_THREADS matrices of different size:
   l, l+delta_l,..., l+(NUM_THREADS-1)delta_l. We create for each one 
   an OpenMP thread. Enabling nested parallelism we can parallelise the 
   diagonalisation of each matrix*/

 if(argc==4)
    {
      sscanf(argv[1],"%d", &l);
      sscanf(argv[2],"%f", &nu);
      sscanf(argv[3],"%d", &NUM_THREADS);
    }

 else
    {
      printf("Some parameters are left...\n");
      exit(0);
    }
 
 /*Migration between cores: disabled*/
 setenv("OMP_PROC_BIND", "TRUE", 1); 
    
 /*Dynamic change of MKL threads: disabled*/
 mkl_set_dynamic(0); 

 /*Nested parallelism: enabled*/
 omp_set_nested(1);

 /*Number of threads employed by OPENMP*/
 /*Nested threads*/
 omp_set_num_threads(NUM_THREADS);   

 /*Number of threads employed by MKL*/
 mkl_set_num_threads(8);

 /* We limit the maximum allowed number of 
    nested, active parallel regions*/
 omp_set_max_active_levels(NUM_THREADS);
    
 delta_l=2000;
  
 #pragma omp parallel for private(arg)
  for(i=0; i<NUM_THREADS; i++)
    { 
     arg=l+i*delta_l; 
     diagonalize(arg);
    }
}


void diagonalize(int l)
{
 int i;
 double *eigen=malloc(l*sizeof(double));
 dcomplex* matrix=malloc((size_t)l*l*sizeof(dcomplex));

 printf("Number of threads employed by MKL: %d\n", mkl_get_max_threads());
  
 printf("%d %1.2f\n", l, nu);
     
 if(matrix == NULL)
  {
   perror("Failed to allocate rows");
  }
 else
  {
   load_matrix(matrix, l);

   printf("l=%d: Matrix loaded\n", l);

   compute_eigenvalues(matrix, eigen, l);

   compute_log_det(eigen, l); 

   write_results(l);
  }

 free(eigen);
 free(matrix);
}

  
void load_matrix(dcomplex *matrix, int l) 
 {
  int n, m;
  FILE *out;
  char name_file[256];  
 
  /*f corresponds to a row of the Toeplitz matrix generate in mathematica
    and saved in a textfile*/ 

  dcomplex* f;
  f=(dcomplex*)malloc(l*sizeof(dcomplex)); 
  
  sprintf(name_file, "matrix_g_%1.2f.dat", nu);

  out=fopen(name_file, "rt");

  for(n=0; n<l; n++)
    {
     fscanf(out,"%*d %lf %lf \n", &f[n].re, &f[n].im); 
    }  
  fclose(out);

  /* * in %*d is works as the assignement-supression character*/
   
  /*In order to construct the submatrix we must follow the directives 
    of the Intel MKL library*/
  

  /*We construct the upper triangular part of the block*/ 
  for(n=0; n<l; n++)
    {
     for(m=n; m<l; m++)
       {
        matrix[n*l+m].re=f[m-n].re;
        matrix[n*l+m].im=f[m-n].im;   
       }
    }

  /*We construct the lower triangular part of the block*/
  for(m=0; m<l-1; m++)
    {
     for(n=0; n<m+1; n++)
       {
        matrix[(m+1)*l+n].re=f[m+1-n].re;
        matrix[(m+1)*l+n].im=-f[m+1-n].im;
       }
     }

 free(f);
}

void compute_eigenvalues(dcomplex *submatrix, double *eigen, int l)
{
 int lwork, info;
 dcomplex* work;
 dcomplex wkopt;
 /* rwork dimension should be at least max(1,3*n-2) */
 double *rwork=malloc((3*l-2)*sizeof(double));

 /* Query and allocate the optimal workspace */
        
 lwork = -1;
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

void compute_log_det(double *eigen, int l)
{
  int n;

  log_det=0;  
 
  for(n=0; n<l; n++)
    log_det+=logl(fabsl(eigen[n]));
     
  printf("l=%d log_det=%lf\n", l, log_det);

}

void write_results(int l)
{
  char name_file[256];

  sprintf(name_file, "sublog_%1.2f.dat", nu);
  out=fopen(name_file, "at");
  fprintf(out,"%d %.10lf\n", l, log_det);
  fclose(out);
}  
