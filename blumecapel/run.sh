#!/bin/bash
rm -r outputs
mkdir outputs
gcc -g -o blmcp blumecapel.c -lm -lgslcblas -lgsl
gcc -g -o bcinf bc_inf.c -lm -lgslcblas -lgsl

declare -i nMeas=40
mu=$(bc<<<"0.25")
declare -i L=8
J=$(bc<<<"1.0")
Tmin=$(bc<<<"0.75")
dT=$(bc<<<"0.05")
Tmax=$(bc<<<"$Tmin+$dT*$nMeas")
p=$(bc<<<"0.0325")
declare -i MCS=100000

echo L = $L, J = $J, mu = $mu, Tmin = $Tmin, Tmax=$Tmax, p = $p,nMeas = $nMeas, MCS = $MCS
for ((c = 0;c<$nMeas;c=c+1)); 
do
    T=$(bc<<<"$Tmin+$c*$dT")
    ./blmcp $L $J $mu $T $p $MCS $c;
    echo -------------------$c----------------------
done;

echo Inizio Inferenza 
./bcinf $L $MCS $Tmin $Tmax $dT $nMeas 
