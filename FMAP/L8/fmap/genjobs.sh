if [ "$1" == "" ];then
	np=`lscpu |grep "^CPU(s):"|awk '{print $2}'`
else
	np=$1
fi

nn=0
for j in `seq 90 -6 6`
do
        for T in 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15
        do
		let nn=nn+1
                #echo $T $j
		echo "bash run.sh $T $j &"
		m=$((nn % np))
		if [ $m -eq 0 ];then
			echo wait
		fi
        done
done
echo wait
