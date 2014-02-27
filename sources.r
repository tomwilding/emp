# Source code
source('aux/sim.r')
source('aux/takeEveryOther.r')
source('aux/breakTime.r')
source('aux/rSquareError.r')
source('aux/sumData.r')

source('core/fitOverTimeMulti.r')
source('core/sseSIRMulti.r')
source('core/evalSIRMulti.r')
source('core/fitInRangeParallel.r')
source('core/setSolver.r')
source('core/sir.r')

source('plot/reconstructPlot.r')
source('plot/plotRSq.r')

# Packages
install.packages("deSolve")
install.packages("doMC")
install.packages("foreach")
install.packages("GillespieSSA")

## Create package (functions and descriptions)
# source('sources.r') # To get dependencies
# package.skeleton("epi") # Construct package - Required to add names

## Install from source using
# install.packages("epi", repos=NULL, type="source")	(To install the package)
# require("epi")	(to load the package)