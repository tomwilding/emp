getFineGrain <- function(array, granArray) {
	mapply(function(x,y) ((x - 1) * (1/y)) + 1, array, granArray)
}

# getIndex <- function(x, gran) {
# 	((x - 1) * (1/gran)) + 1
# }