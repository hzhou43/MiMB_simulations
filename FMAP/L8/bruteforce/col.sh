for T in 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15
do
	cat /dev/null >${T}/res.txt
	for j in `seq -w 6 6 90`
	do
                grep -e "#Density" ${T}/out.0.${j} |awk '{print $2,$8}' >>${T}/res.txt 
        done
done
