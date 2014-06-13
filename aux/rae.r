rae <- function(preds, data) {
	se <- sum(abs(preds - data))
	mean <- mean(data)
	sae <- sum(abs(mean - data))
	rae <- se / sae
}