#Dan's Stop Looper


###Setup

Check out the [cmstas/Software](https://github.com/cmstas/Software) repository someplace. Copy the files `dataMCplotMaker.cc`, `dataMCplotMaker.h`, and `PlotMakingTools.h` from there into this directory.

To compile all the code, simply type `make`. The plotmaking scripts mentioned above were designed to compile in ROOT, so you may find you have to make a few small manual fixes to make them work with g++. The latest version of the plotmaking scripts that I can _guarantee_ to work with no errors/warnings is the version from Nov. 5, 2015.

###Running

To run the entire looper, do `./runLooper`. To see a list of command-line options, do `./runLooper help`.

Stacked histograms will be saved in the 'plots' directory. Raw histograms will be stored in `plots.root`. Datacards will be saved in the 'datacards' directory. The LaTeX code for the yield table(s) and dilepton background estimate table will be printed to the screen.