contourPlot <- function(optimParams, initConds, dataIn, timeIn, confidenceLevels) {
	require(emdbook)

	# Evaluate the negative log likelihood function at range of points
	npoints <- 40
	beta <- optimParams[1]
	gamma <- optimParams[2]
	# Region around optimal parameters to explore
	delta <- 0.5;
	db <- -beta*delta;
	dg <- -gamma*delta;
	# Use relative error in lower and upper bounds
	# b <- seq(from=beta-db, to=beta+db, length=npoints)
	# g <- seq(from=gamma-dg, to=gamma+dg, length=npoints)
	b <- seq(from=-10.5, to=-9.7, length=npoints)
	g <- seq(from=-1.7, to=-0.8, length=npoints)
	# # Beta Uncertainty
	# # z <- lapply(x, sirNegLL, g0=optimParams[2], s0=optimParams[3], initConds=initConds, dataIn=truncData, timeIn=truncTimes)
	# # plot(x,z)

	# Beta and Gamma uncertainty
	z <- apply2d(sirNegLL, b, g, s0=optimParams[3], initConds=initConds, dataIn=dataIn, timeIn=timeIn)
	minz <- min(z)
	# Calculate levels to contour
	percentiles <- sapply(confidenceLevels, function(x) minz+0.5*qchisq(x,2))
	# Plot contour of values for points
	setEPS()
	postscript("graphs/CDCMLE/Contour.eps")	
	image(b,g,z)
	contour(b,g,z, levels=percentiles, add=TRUE)
	# Add contour lines to plot
	points(coef(mle)[1],coef(mle)[2],pch=16)
	title(main="Contour Plot of showing Confidence Levels for CDC Data", cex.main=1, cex.axis=0.8)
	# persp(b,g,z,theta=30, phi=30, expand=0.6,
	# 		col='lightblue', shade=0.75, ltheta=120,
	# 		ticktype='detailed')
	dev.off()
}