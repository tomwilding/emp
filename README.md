Source the R script 
`source("setPackages.r")` 
to install all required packages.

to run the fitting process (For example for processing the Synthetic example)
`source("launch/process<X>.r")`
The specified R file loads the required data source and calls package methods. Where the "<X>" in "launch/process<X>.r" is replaced by the file name for example "launch/processMix.r".

The results of the fitting are stored in the "output/data" directory.

to plot the results of the fitting process (For example for processing the Synthetic example)
`source("plot/plotMix.r")`
The specified R file loads the results of the fitting stage and plots all graphs.

The output of the plots are stored in the "output/graphs" directory (The directory is specified in the "plotConfig" object within the "launch/process<X>.r" file, ensure that the empty directory exists before plotting).

The project directory structure is as follows:
    "launch" - launch fitting process
    "output" - output files and data
    "data" - input data sources
    "plot" - plotting the results

The package "epi" contains the core fitting functions.