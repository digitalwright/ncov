map <- function(tree, bases, ancestor, showPlot=FALSE) {
	# map bases onto the tree
	reps <- 0L
	f <- function(dend) {
		if (is.leaf(dend)) {
			attr(dend, "state") <- bases[dend]
			if (showPlot) {
				if (bases[dend] != ancestor &&
					bases[dend] %in% DNA_BASES) {
					attr(dend, "edgePar") <- list(col="red")
				} else if (bases[dend] == ancestor) {
					attr(dend, "edgePar") <- list(col="blue")
				}
			}
		} else {
			dend[[1]] <- f(dend[[1]])
			dend[[2]] <- f(dend[[2]])
			s1 <- attr(dend[[1]], "state")
			s2 <- attr(dend[[2]], "state")
			if (all(s1 %in% DNA_BASES)) {
				if (all(s2 %in% DNA_BASES)) { # both
					if (any(s1 %in% s2)) { # matching
						if (length(s1) > 1 || length(s2) > 1) {
							t <- table(c(s1, s2))
							attr(dend, "state") <- names(t)[which(t==max(t))]
							
							# coalescence to ancenstral state
							if (length(attr(dend, "state"))==1 &&
								attr(dend, "state")==ancestor) {
								if (showPlot)
									attr(dend, "edgePar") <- list(col="green", lwd=2)
								reps <<- reps + 1L
							}
						} else { # matching
							attr(dend, "state") <- s1
						}
					} else { # difference
						attr(dend, "state") <- unique(c(s1, s2))
					}
				} else { # inherit
					attr(dend, "state") <- s1
				}
			} else {
				if (all(s2 %in% DNA_BASES)) { # inherit
					attr(dend, "state") <- s2
				} else { # neither
					attr(dend, "state") <- "" # unknown
				}
			}
		}
		dend
	}
	new <- f(tree)
	if (showPlot)
		plot(new, leaf="none")
	s1 <- attr(new[[1]], "state")
	s2 <- attr(new[[2]], "state")
	if (all(s1 %in% DNA_BASES) &&
		all(s2 %in% DNA_BASES) &&
		any(!(s1 %in% s2)))
		reps <- reps + 1L
	reps
}
