#' @export
spcfactorforecast <- function(preddata, xdata, startdate, date, h, verbose,
							 p=6, k=8, m=1, factorstartdate,
                             chooseK=c("fixed", "BIC", "IC1", "IC2", "IC3"),
							 center="mean", scale="sd", cachedir="./rescache",
                             screeneddata=NULL)
{
    .genfactorforecast(preddata, xdata, startdate, date, h, verbose,
                        p, k, m, factorstartdate, chooseK, center, scale, cachedir,
                        spcfactors, "spcfactors", screeneddata=screeneddata)
}

#' @export
postspcfactorforecast <- function(preddata, xdata, startdate, date, h, verbose,
							 p=6, k=8, m=1, factorstartdate,
                             chooseK=c("fixed", "BIC", "IC1", "IC2", "IC3"),
							 center="mean", scale="sd", cachedir="./rescache",
                             screeneddata=NULL)
{
    .genfactorforecast(preddata, xdata, startdate, date, h, verbose,
                        p, k, m, factorstartdate, chooseK, center, scale, cachedir,
                        postspcfactors, "postspcfactors", screeneddata=screeneddata)
}

#' @export
spcfactors <- function(X, r, center=mean, scale=NULL, maxit=10000, 
							   convEps = 1e-8, verbose=FALSE,
                               psis=c(seq(0,5,0.1), seq(5.5,10,0.5), 11:20))
{
	.SPCBIC <- function(Z, fhat, lambdahat)
	{
		N <- ncol(Z)
		T <- nrow(Z)
		log((1/(N*T))*sum((Z - fhat%*%t(lambdahat))^2)) + sum(lambdahat!=0)*log(N*T)/(N*T)
	}

	if (!is.null(center)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(center, list(x))))
	if (!is.null(scale)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(scale, list(x))), FUN="/")

	Lambdas <- list()
	Fs <- list()
	BICs <- c()

	for (g in 1:length(psis))
	{
		tmp <- tryCatch(.EstimateSPC(X, r, center=NULL, scale=NULL, 
									 psi=psis[g], maxit=maxit, 
									 convEps = convEps, verbose=verbose), 
						error=function(e) NULL)
		if (is.null(tmp))
		{
			Lambdas[[g]] <- NULL		
			Fs[[g]] <- NULL
			BICs <- c(BICs, Inf)
		}
		else
		{	
			Lambdas[[g]] <- tmp$Lambda		
			Fs[[g]] <- tmp$F
			BICs <- c(BICs, .SPCBIC(X, tmp$F, tmp$Lambda))
		}
	}
	# Determine best estimates
	idx <- which.min(BICs)
	F <- Fs[[idx]]
	Lambda <- Lambdas[[idx]]
	Psi <- psis[idx]
 
	return(list(Lambda=Lambda, F=F, Psi=Psi))
}

#' @export
postspcfactors <- function(X, r, center=mean, scale=NULL, penF=NULL, 
								   penLambda=NULL, 
                                   psis=c(seq(0,5,0.1), seq(5.5,10,0.5), 11:20),
								   maxit=10000, convEps = 1e-8, verbose=FALSE)
{
	if (is.null(penF) || is.null(penLambda))
	{
		tmp <- spcfactors(X, r, center, scale, psis=psis, verbose=verbose)
		penF <- tmp$F
		penLambda <- tmp$Lambda
	}

	if (!is.null(center)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(center, list(x))))
	if (!is.null(scale)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(scale, list(x))), FUN="/")

	unscaledX <- X

	T <- nrow(X)
	N <- ncol(X)

	if (!is.list(penF))
	{
		penF <- list(penF)
		penLambda <- list(penLambda)
	}
	postF <- list()
	postLambda <- list()

	for (jj in 1:length(penF))
	{
		Lambda <- c()
		F <- c()

		for (j in 1:r)
		{
			# Init
			f <- penF[[jj]][,j]
			penl <- penLambda[[jj]][,j] -> l
	
			for (i in 1:maxit)
			{
				# Save old values
				f.old <- f; l.old <- l
		
				# Update
				l <- t(lm(X~0+f)$coef)
				l[penl == 0] <- 0
				f <- X%*%l
				f <- f/sqrt(sum(f^2))
				
				# Check for convergence
				fconv <- sum((f-f.old)^2)
				lconv <- sum((l-l.old)^2)
				if (verbose) cat("Comp", j, "\tF conv:", fconv, "\tL conv", lconv, "\n")
				if (is.na(fconv) || is.na(lconv) || is.infinite(fconv) || is.infinite(lconv)) return(NULL)
				if (fconv < convEps && lconv < convEps) break
				if (i == maxit) stop("MaxIter reached!")
			}
	
			# Scale final estimate
			l <- l/sqrt(sum(l^2))
			f <- unscaledX%*%l
	
			# Save estimates and compute residuals
			F <- cbind(F, f)
			Lambda <- cbind(Lambda, l)
			X <- unscaledX - f%*%t(l)
			unscaledX <- X
		}
		postF[[jj]] <- F
		postLambda[[jj]] <- Lambda
	}
	return(list(Lambda=postLambda[[1]], F=postF[[1]]))
}

.EstimateSPC <- function(X, r, center=mean, scale=NULL, psi=0.1, maxit=10000, convEps = 1e-8, verbose=FALSE)
{
	if (!is.null(center)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(center, list(x))))
	if (!is.null(scale)) X <- sweep(X, 2, apply(X, 2, function(x) do.call(scale, list(x))), FUN="/")

	unscaledX <- X

	Lambda <- c()
	F <- c()

	if (length(psi) == 1) psi <- rep(psi,r)
	else if (length(psi) != r) stop ("Wrong length psi")

	for (j in 1:r)
	{
		# Init
		tmp <- svd(X)
		f <- tmp$u[,1]
		l <- tmp$v[,1]*tmp$d[1]

		for ( i in 1:maxit)
		{
			# Save old values
			f.old <- f; l.old <- l
	
			# Update
			l <- t(X)%*%f
			l <- sign(l)*pmax(abs(l)-psi[j],0)
			f <- X%*%l
			f <- f/sqrt(sum(f^2))
			
			# Check for convergence
			fconv <- sum((f-f.old)^2)
			lconv <- sum((l-l.old)^2)
			if (verbose) cat("Comp", j, "\tF conv:", fconv, "\tL conv", lconv, "\n")
			if (is.na(fconv) || is.na(lconv) || is.infinite(fconv) || is.infinite(lconv)) return(NULL)
			if (fconv < convEps && lconv < convEps) break
			if (i == maxit) stop("MaxIter reached!")
		}

		# Scale final estimate
		l <- l/sqrt(sum(l^2))
		f <- unscaledX%*%l

		# Save estimates and compute residuals
		F <- cbind(F, f)
		Lambda <- cbind(Lambda, l)
		X <- unscaledX - f%*%t(l)
		unscaledX <- X
	}
	return(list(Lambda=Lambda, F=F))
}
