#!/bin/bash
name=POWER
echo $name
for ((j=0; j<=8; j=j+4))
do
	#jdec=$(bc <<<"scale=2; $j/20");
	/anaconda3/bin/python3 social_consensus_structure.py 0 0 $name $j
#	for ((i=0; i<=200; i=i+1))
#	do
#		#for ((k=0; k<=200; k=k+1))
#		#do
#			/anaconda3/bin/python3 social_consensus_structure.py 1 $i $name 0
#		#	echo -n "$k " >> results_quenched_$name.txt
#			echo -n "$i " >> results_quenched_$name.txt
#			echo -n "$j " >> results_quenched_$name.txt
#			/anaconda3/bin/python3 solver_consensus.py -n network_structure_$name.json -i influence_distribution_$name.json -s 200 >> results_quenched_$name.txt

	    	#/anaconda3/bin/python3 solver_consensus.py -mm annealed -n network_structure_$name.json -i influence_distribution_$name.json >> results_annealed_$name.txt
#			if [[ $(( i % 1 )) == 0 ]]
#			then
#				echo "propL = $i, skew = $j" 
#			fi
		#done	
#	done
done
