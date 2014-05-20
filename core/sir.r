sir <- function(time, initConds, params) {
	s <- initConds[1]
	i <- initConds[2]
	r <- initConds[3]
	beta <- logisticTransform(1e-6, params[1], 1)
	gamma <- logisticTransform(1e-4, params[2], 1)
	dS <- -beta*s*i
	dI <- beta*s*i - gamma*i
	dR <- gamma*i
	dt <- 1
	list(c(dS,dI,dR,dt))
	# list(c(dS,dI,dR))
}