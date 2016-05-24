# A library of helper functions

# One level BIC wrapper
.BICwrapper1 <- function(func, name, vals)
{
	bestbic <- Inf
	bestres <- NA
	for (v in vals) {
		l <- list(v)
		names(l) <- name
		res <- do.call(func, l)
		if (res$bic < bestbic) {
			bestbic <- res$bic
			bestres <- res
		}
	}
	return(bestres)
}

# Three level BIC wrapper
.BICwrapper3 <- function(func, names, vals1, vals2, vals3)
{
	bestbic <- Inf
	bestres <- NA
	for (v1 in vals1) for (v2 in vals2) for (v3 in vals3) {
		l <- list(v1, v2, v3)
		names(l) <- names
		res <- do.call(func, l)
		if (res$bic < bestbic) {
			bestbic <- res$bic
			bestres <- res
		}
	}
	return(bestres)
}

#' @export
predmse <- function(res)
{
	forecasts <- truevals <- c()
	for (i in 1:length(res)) {
		forecasts <- rbind(forecasts, res[[i]]$forecasts)
		truevals <- rbind(truevals, res[[i]]$truevals)
	}
	colMeans((forecasts - truevals)^2)
}

#' @export
relmse <- function(res1, res2)
{
	predmse(res1)/predmse(res2)
}

.cachefilename <- function(prefix, ...)
{
	l <- list(...)
	name <- c()
	for (i in 1:length(l)) {
		k <- l[[i]]
		for (j in 1:length(k)) name <- paste(name, k[j], sep="_")
	}
	name <- paste(prefix, name, ".rds", sep="")
	return(name)
}

.getfromcache <- function(filename, dir, func)
{
	if (substr(dir, nchar(dir), nchar(dir)) != "/") dir <- paste(dir, "/", sep="")
	if (!file.exists(dir)) dir.create(dir)
	filename <- paste(dir, filename, sep="")
	if (file.exists(filename)) {
		return(readRDS(filename))
	} else {
		if (!is.null(func)) {
            res <- func()
		    saveRDS(res, file=filename)
        } else {
            res <- NULL
        }
		return(res)
	}
}
