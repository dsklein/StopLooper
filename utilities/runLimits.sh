#!/bin/bash

# Run limits for all the corridor-region datacards
for cardfile in `ls combined_corridor/datacard_corridor_*.txt`
do
	echo Running on card ${cardfile}
	massname=${cardfile##*corridor_}
	combine -M Asymptotic --noFitAsimov -n corrLim ${cardfile} > limits_corridor/limits_corridor_${massname}
done

echo
echo

# Run limits for all the Mlb-binned datacards
for cardfile in `ls combined_baseline/datacard_baseline_*.txt`
do
	echo Running on card ${cardfile}
	massname=${cardfile##*baseline_}
	combine -M Asymptotic --noFitAsimov -n baseLim ${cardfile} > limits_baseline/limits_baseline_${massname}
done

# Now grep those limit files, and grab all the expected limits into one places
grep "50.0%" limits_corridor/limits_corridor*.txt > limits_corridor.txt
grep "50.0%" limits_baseline/limits_baseline*.txt > limits_baseline.txt

