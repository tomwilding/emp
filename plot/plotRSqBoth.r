setEPS()
postscript("output/graphs/blur/rSquarePastBoth.eps", height=3)
rsBlur <- read.csv("rsBlur.csv")
rsBlur2 <- read.csv("rsBlurOT.csv")
rsBlur <- rsBlur[,2]
rsBlur2 <- rsBlur2[,2]
print(length(rsBlur))
print(length(rsBlur2))
plot(1:length(rsBlur2), rsBlur2, type="l", col="steelblue", xlab="Time (Days)", ylab="RSquare", ylim=c(-0.2,1))	
lines(1:length(rsBlur), rsBlur, col="black")
# title("Single and Synthedemic past RSquare for Blurred Lines")
legend("bottomright",c("Synthedemic", "Single"), col=c("black", "steelblue"), cex=0.8, lty=1)
title("Time Optimisation and Parallel Time Search RSquare variation")
dev.off()
