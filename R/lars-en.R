#' @importFrom lars lars
#' @importFrom zoo coredata
#' @export
larsenforecast <- function(preddata, xdata, startdate, date, h, verbose,
                           p=6, k=8, m=1, 
                           chooseK=c("fixed", "BIC", "IC1", "IC2", "IC3"),
                           kappa2=0.25, nvar=30)
{
	forecasts <- ps <- ks <- ms <- c()
	for (j in 1:ncol(preddata))
	{
		# Find the X variable that corresponds to prediction var.
		idx <- colnames(xdata) == colnames(preddata)[j]	
		xvar <- xdata[,idx]
		res <- .larsenforecast(y=preddata[,j], x=xvar, xdata=xdata[,!idx], p=p, k=k,
                              m=m, h=h, startdate=startdate, enddate=date, chooseK=chooseK, kappa2=kappa2, nvar=nvar)
		forecasts <- c(forecasts, res$forecast)
		ps <- c(ps, res$p); ks <- c(ks, res$k); ms <- c(ms, res$m)
	}
	names(forecasts) <- names(ps) <- names(ks) <- names(ms) <- colnames(preddata)
	return(list(forecasts=forecasts, p=ps, k=ks, m=ms))
}

.larsenforecast <- function(y, x, xdata, p, k, m, h, startdate, enddate, chooseK, kappa2, nvar)
{
    # Step 1
    # Partial out lags of y (i.e. x) 
	pfunc <- partial(.arforecast, y=y, x=x, h=h, startdate=startdate, enddate=enddate)
	res <- .BICwrapper1(pfunc, "p", p)

    # Step 2
    # Define inputs for LARS
    X <- window(xdata, start=startdate, end=enddate) 
    X <- head(X, nrow(X)-h)
	n <- ncol(X)
	T <- nrow(X)
	yplus <- c(coredata(res$resid), rep(0, n))
	Xplus <- rbind(scale(coredata(X)), diag(n)*sqrt(kappa2))*((1+kappa2)^(-1/2)) # NOTE SCALE

    # Estimate by LARS
	res <- lars(Xplus, yplus, type="lar", normalize=FALSE)
	coefs <- coef(res)
	idx1 <- rowSums(coefs!=0)==nvar
	idx2 <- coefs[idx1,]!=0
	if (sum(idx2)!=nvar) stop("something is wrong LARS-en")

    # Step 3
    # Estimate factors and forecast
	subX <- window(xdata[,idx2], start=startdate)
	res <- pcfactors(subX, k, "mean", "sd") 
	res2 <- list()
	res2$Factors <- ts(res$F, start=startdate, freq=frequency(xdata))
	res2$Lambda <- res$Lambda

	if (chooseK == "BIC") chosenK <- 1:k else if (chooseK == "fixed") chosenK <- k else stop("not implemented")
	pfunc <- partial(.genfactorforecastinner, y=y, x=x, factors=res2$F, h=h, startdate=startdate, enddate=enddate)
	res <- .BICwrapper3(pfunc, c("p", "k", "m"), p, chosenK, m)
    return(res)
}
