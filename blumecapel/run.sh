#!/bin/bash
gcc -g -o bc1 bc_graph.c -lm -lgslcblas -lgsl
for ((c = 0;c<=40;c=c+1)); 
do
    T=$(bc<<<"0.75+$c*0.05")
    ./bc1 8 1 0.25 $T 0.0325 100000 $c;
done;
