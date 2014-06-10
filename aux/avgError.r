avgError <- function(preds, data){
	avgError <- (sum(preds - data) / length(data))
}