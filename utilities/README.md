#Utilities
A bunch of scripts and macros to help in running the stop analysis

####combine_my_cards.sh
A shell script to combine single-region datacards. Copy this into your `CMSSW/src/HiggsAnalysis/CombinedLimit/` directory, and run it from there: `./combine_my_cards.sh`

####runLimits.sh
A script that loops over all the datacards in two directories (hard-coded) and runs the HiggsCombine tool on each. Also uses grep to aggregate the expected limits. Copy it to the CombinedLimit directory and run from there: `./runLimits.sh`

####parseLimits.py
A python macro that loops over the aggregated limit text files, imports the data, and stores plots in ROOT histograms. `python parseLimits.py`

