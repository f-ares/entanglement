# entanglement
This set of codes is intended for computing the Rényi entanglement 
entropy in a stationary state of a quadratic, homogeneous, fermionic
chain through the two-point correlation matrix. For these states the
correlation matrix is block Toeplitz with symbol a 2x2 matrix.

The entries of a row of this block matrix are calculated with a Mathematica
notebook. The Mathematica files are named starting by "math_*". They can be
ran in the bash command line using the correponding bash script the we
include in each folder of the repository. The Mathematica notebook produces
a plain text file with extension "*.dat" with the entries of a row of the 
correlation matrix. This file can be read with the C programs. These programs
construct the entire correlation matrix from this file, compute 
its spectrum and then the Rényi entanglement entropy, whose value
is written in another "*.dat" file. 

The repository is divided into three folders. The codes in each one are
prepared to compute the entanglement entropy in a different particular 
fermionic chain:

- lrk: Long-Range Kitaev chain with power-law decaying pairings. 
          The C programs compute the spectrum of the matrix using 
          GSL. Those that end with "*_disjoint_number.c" compute 
          the entanglement entropy of "number" disjoint intervals.
          
- moebius: fermionic chains with range of coupling 2. In this case
              the codes are prepared to analyse the Moebius symmetry
              of the entanglement entropy. The C program uses GSL
              to compute the spectrum of the correlation matrix.
              
- sublog: study of the determinant of a Toeplitz matrix with a symbol
             that violates the smoothness condition of the Widom theorem.
             In this case, the C programs employ the Intel MKL library to
             obtain the spectrum. They also employ OpenMP to diagonalise
             several matrix of different length simultaneously. 
            
 Filiberto Ares-Departamento de Física Teórica, Universidad de Zaragoza, 50009 Zaragoza (Spain) ares (at) unizar.es
