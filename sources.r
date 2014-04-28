# Source code
source('aux/findStartOffset.r')
source('aux/simSIR.r')
source('aux/simExp.r')
source('aux/takeEveryOther.r')
source('aux/breakTime.r')
source('aux/rSquareError.r')
source('aux/ssError.r')
source('aux/sumData.r')
source('aux/myMean.r')
source('aux/mySd.r')
source('aux/myMeanDiff.r')
source('aux/mySdDiff.r')
source('aux/logisticTransform.r')
source('aux/logit.r')

source('core/fitOverTimeSingle.r')
source('core/fitOverTimeMulti.r')
source('core/sseMulti.r')
source('core/evalMulti.r')
source('core/fitInRangeParallel.r')
source('core/setSolver.r')
source('core/sir.r')
source('core/expDec.r')

source('plot/reconstructPlot.r')
source('plot/plotRSq.r')
source('plot/plotResiduals.r')
source('plot/plotPred.r')

require(doMC)
require(deSolve)

# Packages
# install.packages("deSolve")
# install.packages("doMC")
# install.packages("foreach")
# install.packages("GillespieSSA")

## Create package (functions and descriptions)
# source('sources.r') # To get dependencies
# package.skeleton("epi") # Construct package - Required to add names

## Install from source using
# install.packages("epi", repos=NULL, type="source")	(To install the package)
# require("epi")	(to load the package)