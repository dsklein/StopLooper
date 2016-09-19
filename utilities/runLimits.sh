#!/bin/bash

# Run limits for all the corridor-region datacards
for cardfile in `ls cards_combined/datacard_corridor_*.txt`
do
	echo Running on card ${cardfile}
	massname=${cardfile##*corridor_}
	combine -M Asymptotic --noFitAsimov -n compLim ${cardfile} > limits_combined/limits_corridor_${massname}
done

echo
echo

# Run limits for all the ICHEP datacards
for cardfile in `ls cards_ichep/datacard_std_T2tt*.txt`
do
	echo Running on card ${cardfile}
	massname=${cardfile##*datacard_std_}
	combine -M Asymptotic --noFitAsimov -n ichepLim ${cardfile} > limits_ichep/limits_ichep_${massname}
done

# Now grep those limit files, and grab all the expected limits into one places
grep "50.0%" limits_combined/limits_corridor*.txt > limits_corridor.txt
grep "50.0%" limits_ichep/limits_ichep*.txt > limits_ichep.txt

