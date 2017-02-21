#Dan's Stop Looper


###Setup

Check out the [cmstas/Software](https://github.com/cmstas/Software) repository someplace. Copy the files `dataMCplotMaker.cc`, `dataMCplotMaker.h`, and `PlotMakingTools.h` from there into this directory.

I'm not the author or maintainer of these plotmaking scripts, so I can't promise that every single version will compile in g++ without errors. The latest version I can guarantee to compile correctly and run without errors is [commit 740e6a6](https://github.com/cmstas/Software/tree/740e6a6da2ebc277cfff81d3b164b485ef80a817/dataMCplotMaker), from February 14, 2017.

You'll probably also need to edit the "includes" in the loopers (`ScanChain.C`, `looperCR2lep.C`, `looperCR0b.C`) to point to a local CORE repository so that they can use the DorkyEventIdentifier and the BadEventFilter.

To compile everything, just do `make`. You can clean up the compiled code with `make clean`.

###Running

To run the entire analysis, do `./runLooper`. To see a list of command-line options, do `./runLooper help`. You can use as many options as you want, e.g. `./runlooper lostlep tables estimate cards`.

If you choose to make them, stacked histograms will be saved in the 'plots' directory, and datacards will be saved in the 'datacards' directory. Raw histograms will be saved in various ROOT files in this directory. The LaTeX code for various tables will be printed to the screen.
