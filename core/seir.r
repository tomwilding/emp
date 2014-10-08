seir <- function(time, initConds, params) {
	s <- initConds[1]
	e <- initConds[2]
	i <- initConds[3]
	r <- initConds[4]
	beta <- exp(params[1])
	gamma <- exp(params[2])
	alpha <- exp(params[3])
	dS <- -beta*s*i
	dE <- beta*s*i - gamma*e
	dI <- gamma*e - alpha*i
	dR <- alpha*i
	# dI <- 0
	# list(c(dS,dI,dR,dI))
	list(c(dS,dE,dI,dR))
}