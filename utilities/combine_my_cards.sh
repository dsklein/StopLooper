#!/bin/bash

# Script to combine datacards of the same SUSY masses
# For each SUSY mass, makes two combined datacards - one with some baseline signal regions, and one with my corridor regions.
# Run this from inside your CMSSW/src/HiggsAnalysis/CombinedLimit/ directory.
# Right now, the script assumes the existence of four subdirectories: uncombined_baseline, uncombined_corridor, combined_baseline, combined_corridor

for susymass in `ls -1 uncombined_baseline/ uncombined_corridor/ | grep -o 'T2tt_[[:digit:]]\+_[[:digit:]]\+' | sort | uniq`
do
	argstring_base=''
	argstring_corr=''
	for filename in `ls uncombined_baseline/datacard_*mlb*_${susymass}.txt`;do argstring_base+=" ${filename}";done
	for filename in `ls uncombined_corridor/datacard_corr*50combo_${susymass}.txt`;do argstring_corr+=" ${filename}";done
	echo Combining cards for ${susymass}...
	python scripts/combineCards.py ${argstring_base} > combined_baseline/datacard_baseline_${susymass}.txt
	python scripts/combineCards.py ${argstring_corr} > combined_corridor/datacard_corridor_${susymass}.txt
done
