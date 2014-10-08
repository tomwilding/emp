wData <- read.csv("data/whitney.csv", header = TRUE)
plotConfig <- list(title="Whitney Huston Downloads", fileName="output/graphs/", dataFile="output/data/blur2/blurData90I.RData", envFile="output/data/blur2/blurEnv.RData", pat=45, rat=300)

times <- wData[,1]
data <- wData[,2]
setEPS()
r <- plotConfig$run
graphName <- "whitney.eps"
postscript(paste(plotConfig$fileName, graphName, sep=''))	
par(mar=c(7.1,4.1,4.1,2.1))
plot(times, data, xlab='Days after death', ylab='Infected Individuals', col='steelblue')
lines(times, data, col='steelblue')
title(main=plotConfig$title, cex.main=1, cex.axis=0.8)
dev.off()