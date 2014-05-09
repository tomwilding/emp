simExp <- function(g,S0){
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

	out <- ssa(x0=c(X=S0), a=c("c*X"), nu=matrix(-1), parms=c(c=g), tf=100, censusInterval=1)
}