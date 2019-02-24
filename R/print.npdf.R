##' Print output from a fitted nonparametric discrete frailty model
##'
##' Prints estimates, standard errors, likelihood and model comparison statistics from a fitted nonparametric discrete frailty model
##'
##' @param x A fitted nonparametric discrete frailty model, as returned by \code{\link{npdf_cox}} with \code{estK=FALSE}
##'
##' @param digits Number of significant digits to present, by default \code{max(1, getOption("digits") - 4)}
##'
##' @param ... Further options (currently unused)
##'
##' @export
print.npdf <- function(x, digits = NULL, ...){
    if (is.null(digits))
        digits <- max(1, getOption("digits") - 4)
    if (!is.null(x$call))
        cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n", sep = "")
    with(x, {
        cat(sprintf("\nNonparametric discrete frailty Cox model fit with K=%s latent populations\n\n", K))
        if( K != 1 ){
          est <- c(p, w[-1], beta)
        }else{
          est <- c( beta )
        }
        cat("Estimated parameters and standard errors:\n")
        res <- data.frame(est=est)
        if (!is.null(x$seLouis)) res$seLouis <- seLouis
        if (!is.null(x$seExact)) res$seExact <- seExact
        if (!is.null(x$seNumeric)) res$seNumeric <- seNumeric
        if( K!=1 ){
            rownames( res )[1:(2*K-1)] =  c(paste("p",1:K, sep=""),
                                            paste( paste("w",2:K, sep=""),"/",
                                                  rep("w1",(K-1)), sep="" ))
        }
        fitstr <- sprintf("Log-likelihood %s, AIC %s, BIC %s", llik, BIC, AIC)
        print(res, digits=digits)
        cat("\n")
        cat("Log-likelihood:", format(llik, digits=digits), "\n")
        cat("AIC:", format(AIC, digits=digits), "\n")
        cat("BIC:", format(BIC, digits=digits), "\n")
        cat("Fitted K:", format(K_fitted, digits=digits), "\n")
    })
    invisible(x)
}

##' Print output from a nonparametric discrete frailty modelling procedure with automatic model selection.
##'
##' Prints the model comparison statistics comparing models with different numbers of latent populations, and prints the estimates from the optimal model according to the criterion that was specified when calling \code{\link{npdf_cox}} (by default, BIC criterion).
##'
##' @param x An object returned by \code{\link{npdf_cox}} with \code{estK=TRUE}, containing a list of fitted nonparametric discrete frailty models
##'
##' @param digits Number of significant digits to present, by default \code{max(1, getOption("digits") - 4)}
##'
##' @param ... Further options (currently unused)
##'
##' @export
print.npdflist <- function(x, digits=NULL, ...){
    if (is.null(digits))
        digits <- max(1, getOption("digits") - 4)
    with(x, {
        cat("\nCall:\n", paste(deparse(call), sep = "\n", collapse = "\n"),
            "\n", sep = "")
        cat("\nModel comparison:\n")
        print(comparison)
        cat("Optimal K:\n")
        print(Kopt)
        crit <- if (criterion=="Laird") "Laird criterion" else criterion
        cat("\nBest model according to ", crit, ":", sep="")
        print(models[[Kopt[criterion]]])
        cat("\nTo examine other models, look at `models` component\n")
    })
    invisible(x)
}
