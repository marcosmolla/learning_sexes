# Calcualte Gini Index
gini <- function(x, na.rm = FALSE){
	if (!is.numeric(x)){
		warning("'x' is not numeric; returning NA")
		return(NA)
	}
	if (!na.rm && any(na.ind <- is.na(x)))
		stop("'x' contain NAs")
	if (na.rm)
		x <- x[!na.ind]
	n <- length(x)
	mu <- mean(x)
	N <- n * (n + 1)
	ox <- x[order(x)]
	dsum <- drop(2*crossprod(1:n,  ox))
	(dsum / (mu * N)) - 1
}

### Test
# variances <- c(10^0,10^.1,10^.2,10^.3,10^.4,10^.5,10^.6,10^.7,10^.8,10^.9,10^1,10^1.1,10^1.2,10^1.3,10^1.4,10^1.5,10^1.6,10^1.7,10^1.8,10^1.9,10^2)
# g <- c()
# for(v in 1:length(variances)){
# k <- (4^2)/variances[v]
# r <- rgamma(n = 100000, shape = k, scale=(theta<-(4/k)) )
# g[v] <- gini(r)
# }
