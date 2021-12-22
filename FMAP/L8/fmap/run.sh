#!/bin/bash

if [ "$1" == "" ];then
	T=0.70
else
	T=$1
fi

if [ "$2" == "" ];then
	k=6
else
	k=$2
fi
mkdir -p $T
cd $T

	n=`echo "5 * $k" |bc`
	rho=0.`printf "%02d" $k`
	../../../src/mcljfft -T $T -N $n -rc 3.0 -L 8.0 -nc 20000001 -dr 0.3 -seed 3 -ne 2000000 -fs 10000 -ts 1 -vscl 1.06  >out.$rho
		let k=i+j
cd ..
