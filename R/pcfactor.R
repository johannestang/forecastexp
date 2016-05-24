#' @export
pcfactorforecast <- function(preddata, xdata, startdate, date, h, verbose,
							 p=6, k=8, m=1, factorstartdate,
                             chooseK=c("fixed", "BIC", "IC1", "IC2", "IC3", "CH"),
							 center="mean", scale="sd", cachedir="./rescache",
                             screeneddata=NULL)
{
    .genfactorforecast(preddata, xdata, startdate, date, h, verbose,
                        p, k, m, factorstartdate, chooseK, center, scale, cachedir,
                        screeneddata=screeneddata)
}

.genfactorforecast <- function(preddata, xdata, startdate, date, h, verbose,
							 p=6, k=8, m=1, factorstartdate,
                             chooseK=c("fixed", "BIC", "IC1", "IC2", "IC3", "CH"),
							 center="mean", scale="sd", cachedir="./rescache",
                             factormeth=pcfactors, cacheprefix="pcfactors",
                             screeneddata=NULL)
{
	chooseK <- match.arg(chooseK)
	forecasts <- ps <- ks <- ms <- c()

    # Screening?
    if (is.null(screeneddata)) {
        fdata <- xdata
    } else {
        cacheprefix <- paste(cacheprefix, "-screened", sep="")
        fdata <- screeneddata
    }

	# Estimate factors and cache the results
	cname <- .cachefilename(cacheprefix, k, date, startdate, center, scale)
	factorres <- .getfromcache(cname, cachedir, 
                               function() factormeth(window(fdata, start=factorstartdate), k, center, scale))
	factors <- ts(factorres$F, start=factorstartdate, freq=frequency(fdata))

	# Determine number of factors
	if (chooseK %in% c("IC1", "IC2", "IC3")) {
		chosenK <- minIC(window(xdata, start=factorstartdate), factorres$F, 
						 factorres$Lambda, chooseK, center, scale)
	} 
	if (chooseK == "fixed") chosenK <- k 
	if (chooseK == "BIC") chosenK <- 1:k
    if (chooseK == "CH") {
    	cname <- .cachefilename(cacheprefix, "CH", k, date, startdate, center, scale)
	    chosenK <- .getfromcache(cname, cachedir, 
                                 function() CHsel(window(xdata, start=factorstartdate), k, center, scale,
                                                  factorres$F, factorres$Lambda)$k)
    }

	for (j in 1:ncol(preddata))
	{
		# Find the X variable that corresponds to prediction var.
		idx <- colnames(xdata) == colnames(preddata)[j]	
		xvar <- xdata[,idx]

		# Partial out parameters
		pfunc <- partial(.genfactorforecastinner, y=preddata[,j], x=xvar, factors=factors, h=h, startdate=startdate, enddate=date)
		
	   	res <- .BICwrapper3(pfunc, c("p", "k", "m"), p, chosenK, m)
		forecasts <- c(forecasts, res$forecast)
		ps <- c(ps, res$p); ks <- c(ks, res$k); ms <- c(ms, res$m)
	}
	names(forecasts) <- names(ps) <- names(ks) <- names(ms) <- colnames(preddata)
	return(list(forecasts=forecasts, p=ps, k=ks, m=ms))
}

.genfactorforecastinner <- function(y, x, factors, p, k, m, h, startdate, enddate)
{
	# Simple sanity check
	if (m < 1 || is.null(factors)) stop("No factors, you shouldn't be using this function!")

    # Do we have k factors?
    if (NCOL(factors) < k) k <- NCOL(factors)
    warning("We do not have k factors!")

	# Build design matrix
	Z <- y 
	if (p > 0) for (i in 1:p) Z <- cbind(Z, lag(x, -(i-1)))
	for (i in 1:m) Z <- cbind(Z, lag(factors[,1:k], -(i-1))) 
	Z <- window(Z, start=startdate, end=enddate)
    colnames(Z) <- c("Y", paste("X", 1:(ncol(Z)-1), sep=""))

	# Create data frame and remove last h rows where y is NA
	Zc <- as.data.frame(Z)
	Zc <- head(Zc, nrow(Zc)-h)

	# Is everything ok? Due to lags of x and the factors we may
    # have no more than max(p-1, m-1) rows with NAs in the beginning.
	if (any(is.na(Zc[max(p, m):nrow(Zc),]))) stop("NAs encountered")
    Zc <- na.omit(Zc)

	# Compute forecast
	model <- lm(as.data.frame(Zc), singular.ok=FALSE)
	frc <- tail(predict(model, newdata=Z), 1)
	bic <- BIC(model)
	return(list(forecast=frc, bic=bic, p=p, k=k, m=m))
}

##########################

#' @export
pcfactors <- function(X, r, center="mean", scale="sd")
{
	if (!is.null(center)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(center, list(x))))
	if (!is.null(scale)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(scale, list(x))), FUN="/")
	tmp <- eigen(t(X)%*%X)
	Lambda <- tmp$vec[,1:r]
	F <- X%*%tmp$vec[,1:r]
	return(list(Lambda=Lambda, F=F))
}

#' @export
minIC <- function(X, F, Lambda, meth=c("IC1", "IC2", "IC3"), center=NULL, scale=NULL)
{
	meth <- paste(".", match.arg(meth), sep="")
	if (!is.null(center)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(center, list(x))))
	if (!is.null(scale)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(scale, list(x))), FUN="/")
    res <- c()
    for (i in 1:ncol(F)) res <- c(res, do.call(meth, list(X,F,Lambda,i)))
	return(which.min(res))
}

.IC1 <- function(X, F, Lambda, k)
{
    N <- ncol(X)
    T <- nrow(X)
    pen <- k*(N+T)/(N*T)*log(N*T/(N+T))
    return(log(sum((X - F[,1:k,drop=FALSE]%*%t(Lambda[,1:k,drop=FALSE]))^2)) + pen)
}

.IC2 <- function(X, F, Lambda, k)
{
    N <- ncol(X)
    T <- nrow(X)
    pen <- k*(N+T)/(N*T)*log(min(N,T))
    return(log(sum((X - F[,1:k,drop=FALSE]%*%t(Lambda[,1:k,drop=FALSE]))^2)) + pen)
}

.IC3 <- function(X, F, Lambda, k)
{
    N <- ncol(X)
    T <- nrow(X)
    pen <- k*(log(min(N,T))/min(N,T))
    return(log(sum((X - F[,1:k,drop=FALSE]%*%t(Lambda[,1:k,drop=FALSE]))^2)) + pen)
}

#' @export
CHsel <- function(X, rmax, center="mean", scale="sd", F=NULL, L=NULL, lambda=NULL, Verbose=FALSE)
{
	T <- nrow(X)
	n <- ncol(X)
	if (!is.null(center)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(center, list(x))))
	if (!is.null(scale)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(scale, list(x))), FUN="/")

	if (is.null(lambda)) lambda <- c(seq(1.5, 2.5, 0.25), 3:15)*min(n,T)^0.45

	nlam <- length(lambda)
	P1 <- ((n+T)/(n*T))*log(n*T/(n+T))
	maxiter <- 5000
	Del <- 2e-3
	tol <- 5e-4
	eta <- 1e-4
	stepscale <- 1
	stepd <- 750

	alfa <- 0.25

	if (is.null(F) | is.null(L))
	{
		tmp <- eigen(X%*%t(X)/(n*T))
		f <- tmp$vec[,1:rmax,drop=FALSE]
		d <- diag(tmp$val[1:rmax])
		Fhat <- f*sqrt(T)
		LhatOLS <- t(X)%*%Fhat/T
	}
	else
	{	
		Fhat <- F
		LhatOLS <- L
	}

	Lhat_new <- list(nlam)
	SSE <- matrix(0, nlam, 1)
	BIC1n <- matrix(0, nlam, 1)
	nf <- matrix(0, nlam, 1)

	for (k in 1:nlam)
	{
		if (Verbose) cat("\nk", k, "\n")
		btahat <- matrix(0, n*rmax, maxiter)
		btahat0 <- t(LhatOLS)
		btahat[,1] <- btahat0
		diff <- 100
		g <- matrix(0, n*rmax, maxiter)
		gn <- matrix(0, n*rmax, maxiter)

		for (i in 1:maxiter)
		{
			if (Verbose) cat("\riter", i)
			# Step 1
			Ltrans <- matrix(btahat[,i], rmax, n)
			vecdiag <- diag(Ltrans%*%t(Ltrans)/n)
			g1 <- t(Fhat)%*%(X-Fhat%*%Ltrans)/n
			g2 <- -lambda[k]/n*Ltrans
			for (j in 1:rmax) g2[j,] <- g2[j,]/(vecdiag[j])^(1-alfa)+eta
			g2 <- alfa*g2

			# Step 2
			g[,i] <- (g1+g2)*(abs(btahat[,i])>tol)

			# Step 3
			gn[,i] <- g[,i]
			if (max(abs(g[,i])) != 0)
			{
				g[,i] <- stepscale*g[,i]/(max(abs(g[,i]))*(i/stepd+1))
			}
			else
			{
				g[,i] <- 0
			}

			# Step 4
			btahat[,i+1] <- btahat[,i] + Del*g[,i]
			diff <- max(abs(btahat[,i+1] - btahat[,i]))

			# Check convergence
			if (diff <= tol) break
			if (i == maxiter) warning("maxiter reached")
		}
		Lhat_new[[k]] <- t(matrix(btahat[,i], rmax, n))
		nf[k] <- sum(colSums(abs(Lhat_new[[k]])>tol)>0)
		SSE[k] <- sum(diag(t(X-Fhat%*%t(Lhat_new[[k]]))%*%(X-Fhat%*%t(Lhat_new[[k]]))))/(n*T)
		BIC1n[k] <- log(SSE[k]) + nf[k]*P1*log(log(n))
	}
	k_opt1n <- which.min(BIC1n)
	if (Verbose) print(nf)
	return(list(k=nf[k_opt1n], F=Fhat, L=Lhat_new[[k_opt1n]]))
}
