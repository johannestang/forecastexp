#' @importFrom glmnet glmnet
#' @export
lassoforecast <- function(preddata, xdata, startdate, date, h, verbose, p=6)
{
    forecasts <- ps <- c()
	for (j in 1:ncol(preddata))
	{
		# Find the X variable that corresponds to prediction var.
		idx <- colnames(xdata) == colnames(preddata)[j]	
		xvar <- xdata[,idx]
		pfunc <- partial(.lassoforecast, y=preddata[,j], x=xvar, xdata=xdata[,!idx], h=h, startdate=startdate, enddate=date)
		res <- .BICwrapper1(pfunc, "p", p)
		forecasts <- c(forecasts, res$forecast)
		ps <- c(ps, res$p)
	}
	names(forecasts) <- names(ps) <- colnames(preddata)
	return(list(forecasts=forecasts, p=ps))
}

.lassoforecast <- function(y, x, xdata, p, h, startdate, enddate)
{
	# Build the design matrix
	Z <- y
	if (p > 0) for (i in 1:p) Z <- cbind(Z, lag(x, -(i-1)))
    Z <- cbind(Z, xdata)
	Z <- window(Z, start=startdate, end=enddate)
    colnames(Z) <- c("Y", paste("X", 1:(ncol(Z)-1), sep=""))

	# Create data frame and remove last h rows where y is NA
	Zc <- as.data.frame(Z)
	Zc <- head(Zc, nrow(Zc)-h)
	if (any(is.na(Zc))) stop("NAs encountered")
	
	# Compute forecast
	model <- glmnet(as.matrix(Zc[,-1]), as.matrix(Zc[,1]), standardize=TRUE, 
				    penalty.factor=c(rep(0,p), rep(1,ncol(xdata))))
	bics <- c()
	for (j in 1:ncol(model$beta)) {
		bics <- c(bics, .lassobic(model$a0[j], model$beta[,j], model$df[j]+1, 
                                  as.matrix(Zc[,1]), as.matrix(Zc[,-1])))
	}
	frc <- predict(model, tail(Z, h)[,-1], s=model$lambda[which.min(bics)])[h]
	return(list(forecast=frc, bic=min(bics), p=p))
}

.lassobic <- function(a0, beta, novars, y, x)
{
	coef <- c(a0, beta)
	resid <- y - cbind(1,x)%*%coef
	rss <- t(resid)%*%resid
	n <- length(y)
	bic <- n*log(rss/n) + novars*log(n)
	bic
}
