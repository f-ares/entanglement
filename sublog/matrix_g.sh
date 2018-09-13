#!/bin/bash                                                                                                                                                 

#This script computes a row of the Toeplitz matrix with symbol g_\nu

l_max=$1 #Dimension of the matrix
nu=$2    #Exponent of the symbol

#Each step of the following loop calls the Mathematica notebook
#math_matrix_g that computes in parallel the entries of 
#the row between old and l

l=$3
delta_l=$4
old=$5

while (( $(echo "$l <= $l_max" | bc)==1)) ; do
    echo $l

    cat math_matrix_g |  sed "s%myl%$l%g" | sed "s%myold%$old%g" | sed "s%mynu%$nu%g" | math

    old=$l+1

    l=$(echo "$l + $delta_l" | bc)
done
