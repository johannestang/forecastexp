#' @export
cfpcforecast <- function(preddata, xdata, startdate, date, h, verbose,
                           p=6, k=8, m=1, 
                           chooseK=c("fixed", "BIC", "IC1", "IC2", "IC3"))
{
	chooseK <- match.arg(chooseK)
	forecasts <- ps <- ks <- ms <- c()
	for (j in 1:ncol(preddata))
	{
		# Find the X variable that corresponds to prediction var.
		idx <- colnames(xdata) == colnames(preddata)[j]	
		xvar <- xdata[,idx]
		res <- .cfpcforecast(y=preddata[,j], x=xvar, xdata=xdata, p=p, k=k,
                              m=m, h=h, startdate=startdate, enddate=date, chooseK=chooseK)
		forecasts <- c(forecasts, res$forecast)
		ps <- c(ps, res$p); ks <- c(ks, res$k); ms <- c(ms, res$m)
	}
	names(forecasts) <- names(ps) <- names(ks) <- names(ms) <- colnames(preddata)
	return(list(forecasts=forecasts, p=ps, k=ks, m=ms))
}

.cfpcforecast <- function(y, x, xdata, p, k, m, h, startdate, enddate, chooseK)
{
	X <- window(xdata, start=startdate, end=enddate)
	y <- window(y, start=startdate)
	Xc <- head(X, nrow(X)-h)
    
    # Step 1
    yhat <- c()
    for (i in 1:ncol(xdata))
    {
		xx <- as.matrix(Xc[,i])
        yy <- as.matrix(y)
		model <- lm(yy~xx)
		xx <- as.matrix(X[,i])
		fit <- cbind(1, xx)%*%model$coef
        yhat <- cbind(yhat, fit)
    } 
    res <- pcfactors(yhat, k, center="mean", scale=NULL)
    factors <- ts(res$F, start=startdate, freq=frequency(xdata))

	# Determine number of factors
	if (chooseK %in% c("IC1", "IC2", "IC3")) {
		chosenK <- minIC(yhat, res$F, 
						 res$Lambda, chooseK, center="mean", scale=NULL)
	} 
	if (chooseK == "fixed") chosenK <- k 
	if (chooseK == "BIC") chosenK <- 1:k

	pfunc <- partial(.genfactorforecastinner, y=y, x=x, factors=factors, h=h, startdate=startdate, enddate=enddate)
	res <- .BICwrapper3(pfunc, c("p", "k", "m"), p, chosenK, m)
    return(res)
}

