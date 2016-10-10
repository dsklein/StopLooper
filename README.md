#Dan's Stop Looper


###Setup

Check out the [cmstas/Software](https://github.com/cmstas/Software) repository someplace. Copy the files `dataMCplotMaker.cc`, `dataMCplotMaker.h`, and `PlotMakingTools.h` from there into this directory.

I'm not the author or maintainer of these plotmaking scripts, so I can't promise that every single version will compile in g++ without errors. The latest version I can _guarantee_ to compile correctly is [commit 758ce41](https://github.com/cmstas/Software/tree/758ce412d19d7b482129f6291168878d2a620a04/dataMCplotMaker), from August 16, 2016.

To compile all the code, simply type `make`.

###Running

To run the entire analysis, do `./runLooper`. To see a list of command-line options, do `./runLooper help`. You can use as many options as you want, e.g. `./runlooper syst lostlep tables estimate cards`.

If you choose to make them, stacked histograms will be saved in the 'plots' directory, and datacards will be saved in the 'datacards' directory. Raw histograms will be saved in various ROOT files in this directory. The LaTeX code for various yield tables and background estimate tables will be printed to the screen.