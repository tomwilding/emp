initConds <- c(500, 1, 1, 0)
params <- c(log(0.005), log(0.1), log(0.1), log(250))
fineTimes <- seq(0, 100, by=0.05);
data <- lsoda(y=initConds, times=fineTimes, func=seir, parms=params)

sus <- data[,2]
exp <- data[,3]
infectious <- data[,4]
rec <- data[,5]

times <- data[,1]
plotConfig <- list(title="SIR Infectious Curve", fileName="output/graphs/", dataFile="output/data/blur2/blurData90I.RData", envFile="output/data/blur2/blurEnv.RData", pat=10, rat=300)

setEPS()
r <- plotConfig$run
graphName <- "seirAll.eps"
postscript(paste(plotConfig$fileName, graphName, sep=''))	
par(mar=c(7.1,4.1,4.1,2.1))
plot(times, infectious, xlab='Days', ylab='Infected Individuals', col='blue', type = "l", ylim=c(0,500))
lines(times, sus, col="green")
lines(times, rec, col="red")
lines(times, exp, col="orange")
title(main=plotConfig$title, cex.main=1, cex.axis=0.8)
legend(x=89.5, y=300, c("S(t)", "E(t)", "I(t)", "R(t)"), col=c("green", "orange", "blue", "red"), lty=c(1, 1, 1), lwd=c(1, 1, 1), cex=0.8)
ParamText <- paste(c("Beta =", ", Gamma =", ", Alpha ="), c(exp(params[1]), exp(params[2]), exp(params[3])), collapse='')
initText <- paste(c("S0 =", ", E0 =", ", I0 =", ", R0 ="), c(initConds[1], initConds[2], initConds[3], initConds[4]), collapse='')
mtext(ParamText, 1, at=plotConfig$pat, padj=8, cex=0.7, col="black")
mtext(initText, 1, at=plotConfig$pat, padj=10, cex=0.7, col="black")
dev.off()