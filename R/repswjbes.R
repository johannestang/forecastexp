#' @export
#' @importFrom xtable xtable
#' @importFrom pryr partial
#' @import macrods
replicate_swjbes <- function()
{
    # Load the dataset and define data function
    data(SW2002, package="macrods")
    SW2002mod <- SW2002
	idx <- SW2002mod$unbal == 0
	SW2002mod$rawdata <- SW2002mod$rawdata[,idx]
	SW2002mod$tcodes <- SW2002mod$tcodes[idx]
    datafunc <- partial(getmacrodata, ds=SW2002mod)

    # Dates:
	dates <- cbind(rep(1970:1998, each=12), rep(1:12,(1998-1970+1)))
	startdate <- fsd <- c(1960,1)

    # Define the variables we wish to forecast and the experiment function
	vars <- c("IP", "GMYXPQ", "MSMTQ", "LPNAG", "PUNEW", "GMDC", "PUXF", "PWFSA")
    expfunc <- partial(forecastexp, datafunc=datafunc, startdate=startdate, 
	    			   dates=dates, vars=vars, verbose=TRUE)

    # Run experiment
    filename <- "results.RData"
    if (file.exists(filename)) 
    {
        load(filename) 
    }
    else 
    {
        res6 <- res12 <- res24 <- list()

        h <- 6
        res6$AR <- expfunc(partial(arforecast, p=0:6), h=h)
        res6$DI <- expfunc(partial(pcfactorforecast, p=0, k=12, chooseK="BIC", factorstartdate=fsd), h=h)
        res6$DIAR <- expfunc(partial(pcfactorforecast, p=0:6, k=12, chooseK="BIC", factorstartdate=fsd), h=h)

        h <- 12
        res12$AR <- expfunc(partial(arforecast, p=0:6), h=h)
        res12$DI <- expfunc(partial(pcfactorforecast, p=0, k=12, chooseK="BIC", factorstartdate=fsd), h=h)
        res12$DIAR <- expfunc(partial(pcfactorforecast, p=0:6, k=12, chooseK="BIC", factorstartdate=fsd), h=h)

        h <- 24
        res24$AR <- expfunc(partial(arforecast, p=0:6), h=h)
        res24$DI <- expfunc(partial(pcfactorforecast, p=0, k=12, chooseK="BIC", factorstartdate=fsd), h=h)
        res24$DIAR <- expfunc(partial(pcfactorforecast, p=0:6, k=12, chooseK="BIC", factorstartdate=fsd), h=h)

        save(res6, res12, res24, file=filename)
    }

    # Make tables
	res <- rbind(
                relmse(res12$DI, res12$AR),
                relmse(res12$DIAR, res12$AR),
                sqrt(predmse(res12$AR))*c(1,1,1,1,12,12,12,12),
                relmse(res6$DI, res6$AR),
                relmse(res6$DIAR, res6$AR),
                sqrt(predmse(res6$AR))*c(1,1,1,1,6,6,6,6),
                relmse(res24$DI, res24$AR),
                relmse(res24$DIAR, res24$AR),
                sqrt(predmse(res24$AR))*c(1,1,1,1,24,24,24,24)
                )
	rownames(res) <- c("DI, h=12", "DI-AR, h=12", "RMSE AR, h=12", 
					   "DI, h=6",  "DI-AR, h=6",  "RMSE AR, h=6",
					   "DI, h=24", "DI-AR, h=24", "RMSE AR, h=24")
    print(res)
	print(xtable(res, digits=4), size="tiny", file="results.tex")
}
