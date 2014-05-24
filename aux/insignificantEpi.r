insignificantEpi <- function(params, epiTypes, k) {
	for (i in 1:k) {
		# Get sub epidemic type and parameters
		subEpiNumParamsOffset <- 0
		subEpiNumParams <- epiTypes[i]
		if (subEpiNumParams > 0) {
			paramsMulti <- params[(subEpiNumParamsOffset + 1) : (subEpiNumParamsOffset + subEpiNumParams)]
			initCondsMulti <- initConds[(subEpiNumParamsOffset + 1) : (subEpiNumParamsOffset + subEpiNumParams)]
			subEpiNumParamsOffset <- subEpiNumParamsOffset + subEpiNumParams

			if(subEpiNumParamsOffset == 4) {
				if(paramsMulti[3] < -5) {
					print("Insignificant Epi")
					
				}
			}
}