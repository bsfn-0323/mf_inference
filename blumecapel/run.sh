#!/bin/bash
mkdir outputs
gcc -g -o blmcp blumecapel.c -lm -lgslcblas -lgsl
gcc -g -o bcinf bc_inf.c -lm -lgslcblas -lgsl
for ((c = 0;c<=40;c=c+1)); 
do
    T=$(bc<<<"0.75+$c*0.05")
    ./blmcp 8 1 0.25 $T 0.0325 100000 $c;
    echo -------------------$c----------------------
done;
./bcinf 8 100000 0.75 2.75 0.05 41 0
