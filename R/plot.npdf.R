##' Survival curves from a fitted nonparametric discrete frailty model
##'
##' This function plots Kaplan-Meier estimates of the survival for each fitted latent population. 
##'
##' @param x A fitted nonparametric discrete frailty model, as returned by \code{\link{npdf_cox}} with \code{estK=FALSE}
##'
##' @param survfit_opts List of arguments to pass to \code{\link{survfit.formula}}
##'
##' @param ... Arguments to pass to \code{\link{plot.survfit}}
##'
##' @examples ## TODO
##' 
plot.npdf <- function(x, survfit_opts = NULL, ...){
    with(x,{
        Y <- model.extract(mf, "response")
        Y <- as.data.frame(unclass(Y))
        groups <- model.extract(mf, "groups")
        Y$pop <- rep(belonging, table(groups))
        survfit.call <- list(formula = Surv(time, status) ~ pop,
                             data = Y)
        surv <- do.call("survfit", c(survfit.call, survfit_opts))
        plot(surv, ...)
    })
    invisible()
}


##' Survival curves from a nonparametric discrete frailty model chosen by automatic model selection
##'
##' @param x An object returned by \code{\link{npdf_cox}} with \code{estK=TRUE}, containing a list of fitted nonparametric discrete frailty models
##'
##' @param K The number of latent populations which identifies the specific model to plot survival curves for.  By default this is the model selected by the criterion specified when calling \code{\link{npdf_cox}}.
##'
##' @param ... Arguments to pass to \code{\link{plot.survfit}}
##'
##' @examples ## TODO
##'
##' 
plot.npdflist <- function(x, K=NULL, ...){
    with(x, {
        if (is.null(K))
            mod <- models[[Kopt[criterion]]]
        mod$mf <- mf
        plot(mod, ...)
    })
    invisible()
}
