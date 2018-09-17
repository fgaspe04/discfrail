##' Survival curves from a fitted nonparametric discrete frailty model
##'
##' This function plots Kaplan-Meier estimates of the survival for each fitted latent population.
##'
##' @param x A fitted nonparametric discrete frailty model, as returned by \code{\link{npdf_cox}} with \code{estK=FALSE}
##'
##' @param NA_estimate binary variable. If TRUE group-specific Nelson-Aalen estimators are plotted, if FALSE group-specific Kaplan-Meier estimators are plotted. The deafult option is FALSE.
##'
##' @param survfit_opts List of arguments to pass to \code{\link{survfit.formula}}
##'
##' @param ... Arguments to pass to \code{\link{plot.survfit}}
##'
##' @examples
##' load('weibdata.rda')
##'
##' dim(weibdata)
##' result = npdf_cox( Surv(time, status) ~ x, groups = family, data = weibdata, K = 4, estK = F, eps_conv=10^-4)
##' plot( result, NA_estimate = T, col = result$belonging )
##'
##' @export

plot.npdf <- function(x, NA_estimate = F, survfit_opts = NULL, ...){
    with(x,{
        Y <- model.extract(mf, "response")
        Y <- as.data.frame(unclass(Y))
        groups <- model.extract(mf, "groups")
        if( NA_estimate == T ){
          plot( NULL, xlim = c(0,150), ylim= c(0,1), xlab = 'Time [days]', ylab = 'Nelson-Aalen estimator'  )
          for( i in 1 : length( belonging ) )
          {
            NA_est = survfit( Surv(time, status) ~ 1 , data = Y[ which( groups ==  i ), ], type="fh")
            y0 = c( 0, -log( NA_est$surv ) )
            sfun0 = stepfun( NA_est$time, y0, f = 0 )
            lines( sfun0, col = belonging[ i ], cex.axis = 1.2, cex.lab = 1.2, do.points = F )
          }

        }else{
          survfit.call <- list(formula = Surv(time, status) ~ groups,
                             data = Y)
          surv <- do.call("survfit", c(survfit.call, survfit_opts))
          plot(surv, ...)
        }
    })
    invisible()
}


##' Survival curves from a nonparametric discrete frailty model chosen by automatic model selection
##'
##' @param x An object returned by \code{\link{npdf_cox}} with \code{estK=TRUE}, containing a list of fitted nonparametric discrete frailty models
##'
##' @param K The number of latent populations which identifies the specific model to plot survival curves for.  By default this is the model selected by the criterion specified when calling \code{\link{npdf_cox}}.
##'
##' @param K The number of latent populations which identifies the specific model to plot survival curves for.  By default this is the model selected by the criterion specified when calling \code{\link{npdf_cox}}.
##'
##' @examples ## TODO
##'
##' @export

plot.npdflist <- function(x, K=NULL, ...){
    with(x, {
        if (is.null(K))
            mod <- models[[Kopt[criterion]]]
        mod$mf <- mf
        plot(mod, ...)
    })
    invisible()
}


nelsonaalen.npdf <- function(x){
    Y <- model.extract(x$mf, "response")
    Y <- as.data.frame(unclass(Y))
    groups <- model.extract(x$mf, "groups")
    Y$pop <- rep(x$belonging, table(groups))
    ngrp <- length(unique(groups))
    res <- vector(ngrp, mode="list") 
    for (i in seq(length=ngrp)){
        est <- survfit(Surv(time, status) ~ 1 , data = Y[groups ==  i, ], type="fh")
        est$y0 <- c( 0, -log( est$surv ) )
        est$sfun0 <- stepfun( est$time, est$y0, f = 0 )
        res[[i]] <- est
    }
    attr(res, "belonging") <- x$belonging
    res
}

##'  x <- npdf_cox( Surv(time, status) ~ x, groups=family, data=weibdata, K = 4, estK=FALSE, eps_conv=10^-4)
##' object <- nelsonaalen.npdf(x)
##' plot.nelsonaalen.npdf(object)
##' plot.nelsonaalen.npdf(object, xlim=c(0,200), ylim=c(0,2), xlab="Follow-up days", ylab="Nelson-Aalen estimate", cols=ifelse(x$belonging==2, "purple", "black"))
##' plot.npdf.tmp(x, type="na")
##' plot.npdf.tmp(x, type="na", na_opts = list(xlim=c(0,200), ylim=c(0,2), xlab="Follow-up days", ylab="Nelson-Aalen estimate", cols=ifelse(x$belonging==2, "purple", "black")))
##' plot.npdf.tmp(x, type="km")
##'
##' 

plot.nelsonaalen.npdf <- function(object, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, cols=NULL, ...){
    ngrp <- length(object)
    sfuns <- lapply(object, "[[", "sfun0")
    if (is.null(xlim)) xlim <- range(lapply(sfuns, knots))
    if (is.null(ylim)) ylim <- c(0,1)
    if (is.null(xlab)) xlab <- 'Time'
    if (is.null(ylab)) ylab <- 'Cumulative hazard'
    if (is.null(cols)) cols <- attr(object, "belonging")
    plot( NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
    for (i in seq(length=ngrp)){
        lines( sfuns[[i]], col = cols[ i ], do.points=FALSE )
    }
}

### TODO colour KM by belonging 
### documentation
### separate function to calculate KM

plot.npdf.tmp <- function(x, type="km", survfit_opts = NULL, na_opts = NULL, ...){
    type <- match.arg(type, c("km", "na"))
    with(x,{
        if( type=="na"){
            object <- nelsonaalen.npdf(x)
            pnacall <- list(object) 
            do.call("plot.nelsonaalen.npdf", c(pnacall, na_opts))
        }else{
            Y <- model.extract(mf, "response")
            Y <- as.data.frame(unclass(Y))
            groups <- model.extract(mf, "groups")
            survfit.call <- list(formula = Surv(time, status) ~ groups,
                                 data = Y)
            surv <- do.call("survfit", c(survfit.call, survfit_opts))
            plot(surv, ...)
        }
    })
    invisible()
}

