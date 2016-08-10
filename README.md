#Dan's Stop Looper


###Setup

Check out the [cmstas/Software](https://github.com/cmstas/Software) repository someplace. Copy the files `dataMCplotMaker.cc`, `dataMCplotMaker.h`, and `PlotMakingTools.h` from there into this directory.

To compile all the code, simply type `make`. The plotmaking scripts mentioned above aren't necessarily designed to compile in g++, so you may find you have to edit them by hand to make them compliant.

###Running

To run the entire looper, do `./runLooper`. To see a list of command-line options, do `./runLooper help`.

Stacked histograms will be saved in the 'plots' directory. Raw histograms will be stored in `plots.root`. Datacards will be saved in the 'datacards' directory. The LaTeX code for the yield table(s) and dilepton background estimate table will be printed to the screen.