#Computes the entanglement entropy after a Lorentz transformation of parameter zeta

#!/bin/bash

l=$1 #Length of the interval
delta_l=$2
l_max=$3

zeta=$4
delta_zeta=$5
zeta_max=$6

alpha=$7 #Parameter of the Rényi entropy

while (( $(echo "$zeta <= $zeta_max" | bc)==1)) ; do

   old=0
   
   echo zeta=$zeta
   
   rm moebius_correl.dat
   touch moebius_correl.dat

   while (( $(echo "$l <= $l_max" | bc)==1)) ; do
     
     echo l=$l
   
     cat math_moebius | sed "s%myzeta%$zeta%g" | sed "s%myl%$l%g" | sed "s%old_l%$old%g" | math

     ./moebius $l $zeta $alpha

     old=$l;
  
     l=$(echo "$l + $delta_l" | bc)
   
   done 
   
   l=$1
   zeta=$(echo "$zeta+$delta_zeta" | bc)
   
done
