#!/bin/bash

for (( i = 1000 ; i <= 16384*1000 ; i*=2))
do
	echo " loop $i"
	time ./harmonic_different_n $i
done  
