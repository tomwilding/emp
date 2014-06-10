# Source code
source('aux/findStartOffset.r')
source('aux/simSIR.r')
source('aux/simExp.r')
source('aux/takeEveryOther.r')
source('aux/breakTime.r')
source('aux/rSquareError.r')
source('aux/ssError.r')
source('aux/relSSError.r')
source('aux/rmse.r')
source('aux/avgError.r')
source('aux/sumData.r')
source('aux/myMean.r')
source('aux/mySd.r')
source('aux/myMeanDiff.r')
source('aux/mySdDiff.r')
source('aux/orderOf.r')
source('aux/testParams.r')
source('aux/logistic.r')
source('aux/logit.r')
source('aux/fitARModel.r')

source('core/fitOverTimeSingle.r')
source('core/fitOverTimeMulti.r')
source('core/detectOutbreak.r')
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
source('plot/analysis.r')

# Packages
# install.packages("deSolve")
# install.packages("doMC")
# install.packages("foreach")
# install.packages("GillespieSSA")
# install.packages("bbmle")
# install.packages("forecast")

## Create package (functions and descriptions)
# source('sources.r') # To get dependencies
# package.skeleton("epi") # Construct package - Required to add names

## Install from source using
# install.packages("epi", repos=NULL, type="source")	(To install the package)
# require("epi")	(to load the package)