#!/bin/bash
name=$1
echo $name
python priority_distribution.py 0 $name 0
for ((i=0; i<=20; i=i+1))
do
	python priority_distribution.py 1 $name $i
	python solver_consensus.py  -pf priority_dict_$name.json -i influence_distribution_$name.json >> results_annealed_$name.txt
	if [[ $(( i % 5 )) == 0 ]]
	then
		echo $i
	fi
done
