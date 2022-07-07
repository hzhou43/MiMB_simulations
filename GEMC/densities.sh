#!/bin/bash

rm -rf phase.txt
for i in $(seq 0.65 0.05 0.65); do
	awk 'NR>50000' T${i}/output.txt | awk '{s1+=$8;s2+=$9}END{print '${i}', s1/NR, s2/NR}' >> phase.txt
done
