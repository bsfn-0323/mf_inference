#!/bin/bash
mkdir outputs
gcc -g -o blmcp bc_graph.c -lm -lgslcblas -lgsl
gcc -g -o inferenza inferenza.c -lm -lgslcblas -lgsl
for ((c = 0;c<=40;c=c+1)); 
do
    T=$(bc<<<"0.75+$c*0.05")
    ./blmcp 8 1 0.25 $T 0.0325 100000 $c;
done;
./inferenza 8 100000 0.75 2.75 0.05 41 0
