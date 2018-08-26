# It computes the entanglement entropy for an interval in the ground state of 
# the Long-Range Kitaev chain

#!/bin/bash
  
l=$1 # Initial length of the interval
l_max=$2 # Maximum length of the interval
delta_l=$3 # Increase in the length of the interval 

delta=$4 # Dumping exponent of the pairings
h=$5 # Chemical potential 

alpha=$6 #Parameter of the Rényi entropy

old=0

rm lrk_correl.dat
touch lrk_correl.dat

while (( $(echo "$l <= $l_max" | bc)==1)) ; do
  echo $l
  
  cat math_lrk_correl | sed "s%myl%$l%g" | sed "s%old_l%$old%g" | sed "s%mydelta%$delta%g" | sed "s%myh%$h%g" | math

  ./lrk_diagonal $l $delta $h $alpha

  old=$l;
  
  l=$(echo "$l + $delta_l" | bc)
done
