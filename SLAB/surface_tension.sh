#!/bin/bash

rm -rf gammaKB.txt gammaIK.txt

for i in $(seq 0.65 0.05 1.10); do
	awk 'NR>=50000' T${i}/log-output_equil* | awk '{s1+=$4; s2+=$5; s3+=$6}END{ print '${i}', 75*(s3/NR - 0.5*(s1+s2)/NR)}' >> gammaKB.txt
	awk 'NR>0' T${i}/histograms.txt | awk '{s1+=$5}END{print '${i}', 0.1*s1/2.0}' >> gammaIK.txt
done
