#!/bin/bash
name=$1
echo $name
~/anaconda3/bin/python3 social_consensus_new.py 0 0 $name
for ((i=0; i<=100; i=i+5))
do
	~/anaconda3/bin/python3 social_consensus_new.py 1 $i $name
	~/anaconda3/bin/python3 solver_consensus.py  -n network_structure_$name.json -i influence_distribution_$name.json >> results_$name.txt
	if [[ $(( i % 5 )) == 0 ]]
	then
		echo $i
	fi
done
