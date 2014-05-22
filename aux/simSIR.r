simSIR <- function(b,g,S0,I0){
	require(GillespieSSA)

	# ===============================================================================
	# Kermack-McKendric SIR model (Brown & Rothery, 1993)
	#===============================================================================

	# The Kermack-McKendrick SIR model is defined as
	# dS/dt = -beta*N*S
	# dI/dt = beta*N*S - gamma*I
	# dR/dt = gamma*I
	#
	# This model consists of two reactions with the following per capita rates,
	# transmission: beta
	# recovery:     gamma

	# Define parameters
	parms <- c(beta=b, gamma=g)

	# Define system
	x0 <- c(S=S0, I=I0, R=0)                      # Initial state vector
	nu <- matrix(c(-1,0,1,-1,0,1),nrow=3,byrow=T) # State-change matrix
	a  <- c("beta*S*I", "gamma*I")                # Propensity vector
	tf <- 1000                                     # Final time
	simName <- "sir"

	# Run the simulations
	# nf <- layout(matrix(c(1,2,3,4),ncol=2, byrow=T))

	# Direct method
	set.seed(2)
	out <- ssa(x0,a,nu,parms,tf,method="D",simName,censusInterval=1)
	# plot(out$data[,1], out$data[,3])
	# print(out$data[,3])
	# readline()
}