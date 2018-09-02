# entanglement
This set of codes is intended for computing the Rényi entanglement 
entropy in a stationary state of a quadratic, homogeneous, fermionic
chain through the two-point correlation matrix. For these states the
correlation matrix is block Toeplitz with symbol a 2x2 matrix.

The entries of a row of this block matrix are calculated with a Mathematica
notebook. The Mathematica files are named starting by math_*. They can be
ran in the bash command line using the correponding bash script the we
include in each folder of the repository. The Mathematica notebook produces
a plain text file with extension .dat with the entries of a row of the 
correlation matrix. This file can be read with the C programs. These programs
construct the entire correlation matrix from this file, compute 
its spectrum and then the Rényi entanglement entropy, whose value
is written in another *.dat file. 

 
