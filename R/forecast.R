#' forecastexp: General purpose forecast experiment framework
#'
#' @docType package
#' @name forecastexp-package
#' @importFrom pryr partial
#' @importFrom parallel mclapply
NULL

#' @export
forecastexp <- function(forecastfunc, datafunc, startdate, dates, vars, h, verbose=FALSE, screen=FALSE, mcore=TRUE)
{
	# Do not compute forecasts for the last h dates as we will have no 
	# true values to compare forecast to.
	fdates <- head(dates, nrow(dates)-h)
	
	# Compute forecasts 
    ffunc <- function(i, datafunc, forecastfunc, startdate, fdates, vars, h, verbose) {
			    if (verbose) cat("Iteration", i, "\n")
			    preddata <- datafunc(fdates[i,], h=h, vars=vars)
			    xdata <- datafunc(fdates[i,])
                if (screen) {
			        screeneddata <- datafunc(fdates[i,], screen=TRUE)
			        res <- forecastfunc(preddata, xdata, startdate, fdates[i,], h, verbose, screeneddata=screeneddata)
                } else {
			        res <- forecastfunc(preddata, xdata, startdate, fdates[i,], h, verbose)
                }
			    return(res)	
		    } 
    if (mcore) {
        results <- mclapply(1:nrow(fdates), ffunc, datafunc, forecastfunc, startdate, fdates, vars, h, verbose)
    } else {
        results <- lapply(1:nrow(fdates), ffunc, datafunc, forecastfunc, startdate, fdates, vars, h, verbose)
    }

	# Compute true values for predicted variables and augment results list.
	ndates <- dates[-c(1:h),]
    tvfunc <- function(i, datafunc, ndates, vars, h, results) {
			    preddata <- datafunc(ndates[i,], h=h, vars=vars)
			    tvs <- tail(preddata, 1)[1,]
			    res <- results[[i]]
			    res$truevals <- tvs
			    return(res)
		    } 
    if (mcore) {
	    results <- mclapply(1:nrow(ndates), tvfunc, datafunc, ndates, vars, h, results)
    } else {
	    results <- lapply(1:nrow(ndates), tvfunc, datafunc, ndates, vars, h, results)
    }
	return(results)
}
