#' @export
arforecast <- function(preddata, xdata, startdate, date, h, verbose, p=6)
{
    forecasts <- ps <- c()
	for (j in 1:ncol(preddata))
	{
		# Find the X variable that corresponds to prediction var.
		idx <- colnames(xdata) == colnames(preddata)[j]	
		xvar <- xdata[,idx]
		pfunc <- partial(.arforecast, y=preddata[,j], x=xvar, h=h, startdate=startdate, enddate=date)
		res <- .BICwrapper1(pfunc, "p", p)
		forecasts <- c(forecasts, res$forecast)
		ps <- c(ps, res$p)
	}
	names(forecasts) <- names(ps) <- colnames(preddata)
	return(list(forecasts=forecasts, p=ps))
}

.arforecast <- function(y, x, p, h, startdate, enddate)
{
	# If p is 0 there are no regressors:
	if (p == 0)
	{
		Z <- window(y, start=startdate)
		model <- lm(Z~1, singular.ok=FALSE)
		frc <- tail(predict(model), 1)
		bic <- BIC(model)
        resid <- ts(model$resid, start=time(Z)[1], freq=frequency(y))
	}
	else {
		# Otherwise build the design matrix
		Z <- y
		for (i in 1:p) Z <- cbind(Z, lag(x, -(i-1)))
		Z <- window(Z, start=startdate, end=enddate)
        colnames(Z) <- c("Y", paste("X", 1:(ncol(Z)-1), sep=""))

		# Create data frame and remove last h rows where y is NA
		Zc <- as.data.frame(Z)
		Zc <- head(Zc, nrow(Zc)-h)
		if (any(is.na(Zc))) stop("NAs encountered")
	
		# Compute forecast
		model <- lm(Zc, singular.ok=FALSE)
		frc <- tail(predict(model, newdata=Z), 1)
		bic <- BIC(model)
        # NOT VALID IF WE REMOVE OBS DUE TO LAGS
        resid <- ts(model$resid, start=time(Z)[1], freq=frequency(y))
	}
	return(list(forecast=frc, bic=bic, p=p, resid=resid))
}

