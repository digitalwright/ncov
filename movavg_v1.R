movavg <- function(x, a=0.01) {
	pad <- 1000
	
	w <- which(is.na(x))
	if (length(w) > 0)
		x[w] <- 0
	
	x <- c(rep(0, pad), x, rep(0, pad))
	
	s1 <- s2 <- numeric(length(x))
	
	for (i in seq_along(x)) {
		if (i == 1) {
			s1[i] <- x[i]
		} else {
			s1[i] <- a*x[i] + (1 - a)*s1[i - 1]
		}
	}
	
	x <- rev(x)
	for (i in seq_along(x)) {
		if (i == 1) {
			s2[i] <- x[i]
		} else {
			s2[i] <- a*x[i] + (1 - a)*s2[i - 1]
		}
	}
	s2 <- rev(s2)
	
	s <- (s1 + s2)/2
	s[-c(seq_len(pad), (length(x) - pad + 1):length(x))]
}
