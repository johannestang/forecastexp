#' @export
ladfactorforecast <- function(preddata, xdata, startdate, date, h, verbose,
							 p=6, k=8, m=1, factorstartdate,
                             chooseK=c("fixed", "BIC", "IC1", "IC2", "IC3"),
							 center="median", scale="sd", cachedir="./rescache",
                             screeneddata=NULL)
{
    .genfactorforecast(preddata, xdata, startdate, date, h, verbose,
                        p, k, m, factorstartdate, chooseK, center, scale, cachedir,
                        ladfactors, "ladfactors", screeneddata=screeneddata)
}

#' @importFrom quantreg rq.fit.br
#' @export
ladfactors <- function(X, r, center=median, scale=NULL, eps=1e-12, maxIter=1000, verbose=FALSE)
{
	# Wrapper function for rq. Let's us do multiple regressions for matrix y.
	.rqq<- function(y,x,size)
	{
		res <- numeric(size)*NA
		for (j in 1:size) res[j] <- rq.fit.br(x, y[,j], tau=0.5)$coef
		return(res)
	}
	# Compute smoothed approximative objective function for a single factor with smoothing par. d
	.objfuncS <- function(par, d=1/(N*T))
	{
		f <- par[1:T]; l <- par[(T+1):(T+N)]
		return(mean(sqrt((X - f%*%t(l))^2 + d^2)))
	}
	# Compute gradient for the smoothed approximative objective function for a single factor with smoothing par. d
	.gradientS <- function(par, d=1/(N*T))
	{
		f <- par[1:T]; l <- par[(T+1):(T+N)]
		tmp <- -(X - f%*%t(l))/(sqrt((X - f%*%t(l))^2 + d^2))
		tlambda <- t(l%x%t(rep(1,T)))
		resF <- rowSums(tmp*tlambda)/(N*T)
		tF <- f%x%t(rep(1,N))
		resLambda <- colSums(tmp*tF)/(N*T)
		return(c(resF, resLambda))
	}

	.optimSLAD <- function(l=NULL, f=NULL)
	{
		control = list(maxit=5000000)

		# Starting values:
		if (is.null(l) && is.null(f))
		{
			l <- rep(1,N)
			f <-rep(1,T)
		}
		convergence <- FALSE
		res <- optim(c(f, l), .objfuncS, gr=.gradientS, method="L-BFGS-B", control=control)
		if (res$convergence != 0) cat("\nERROR in estimation: ", res$convergence, 
									  "\nIter: ", res$iterations, "\nEval: ", res$evaluations, 
									  "\nMessage: ", res$message, "\n")
		if (verbose) cat(res$message, "\n")
		if (res$convergence == 0) convergence <- TRUE
		f <- res$par[1:T]
		l <- res$par[(T+1):(T+N)]
		tmp <- sqrt(sum(l^2))
		l <- l/tmp
		f <- f*tmp
		return(list(convergence=convergence, f=f, l=l))
	}
	.nipals <- function(f,s)
	{
		convergence <- FALSE
		error <- FALSE
		for (i in 1:maxIter)
		{
			fold <- f				
			oldVal <- t(f)%*%f
			tryCatch(.rqq(X,f, N) -> l, error=function(e) {error <- TRUE})
			l <- l/sqrt(sum(l^2))
			tryCatch(.rqq(t(X),l, T) -> f, error=function(e) {error <- TRUE})
			if (any(is.na(f)) || any(is.na(l)) || sd(f) < 1e-12 || sd(l) < 1e-12) error <- TRUE
			conv <- abs(t(f)%*%f - oldVal)
			if (conv < eps) convergence <- TRUE
			if (verbose) cat("Component", s, "\tIter", i, "\tConv", conv, "\tFconv", sum((f-fold)^2), "\tT", nrow(X), "\tN", ncol(X), "\n")
			if (error || convergence) break
		}
		return(list(convergence=convergence, f=f, l=l))
	}

	# Center X and determine its dimensions
	X <- as.matrix(X) # Important!! No dataframes please!
	
	if (!is.null(center)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(center, list(x))))
	if (!is.null(scale)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(scale, list(x))), FUN="/")

	T <- nrow(X)
	N <- ncol(X)

	# F and Lambda will hold the final estimates. Let them be empty for now.
	F <- c()
	Lambda <- c()	

	error <- FALSE

	# Loop over the number of components to estimate
	for (s in 1:r)
	{	
		tmp <- .optimSLAD()
		if (!tmp$convergence) error <- TRUE else 
		{ 
			tmp <- .nipals(tmp$f, s)
			if (!tmp$convergence) error <- TRUE else { f <- tmp$f; l <- tmp$l }
		}

		if (error) break else 
		{
			F <- cbind(F, f)
			Lambda <- cbind(Lambda, l)
			if (s < r) X <- X - f%*%t(l)
		}
	}
	return(list(F=F, Lambda=Lambda))
}

