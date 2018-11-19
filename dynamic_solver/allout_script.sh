#!/bin/bash
name=$1
echo $name
~/anaconda3/bin/python3 social_consensus_structure.py 0 0 $name
for ((i=0; i<=20; i=i+5))
do
	~/anaconda3/bin/python3 social_consensus_structure.py 1 $i $name
	~/anaconda3/bin/python3 solver_consensus.py  -n network_structure_$name.json -i influence_distribution_$name.json --allout -o results_quenched_$name\_$i.json
    ~/anaconda3/bin/python3 solver_consensus.py -mm annealed -n network_structure_$name.json -i influence_distribution_$name.json --allout -o results_annealed_$name\_$i.json
	if [[ $(( i % 5 )) == 0 ]]
	then
		echo $i
	fi
done
