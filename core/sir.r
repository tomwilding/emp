sir <- function(time, initConds, params) {
	s <- initConds[1]
	i <- initConds[2]
	r <- initConds[3]
	beta <- exp(params[1])
	gamma <- exp(params[2])
	dS <- -beta*s*i
	dI <- beta*s*i - gamma*i
	dR <- gamma*i
	list(c(dS,dI,dR))
}