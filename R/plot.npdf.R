##' Survival or cumulative hazard curves from a fitted nonparametric discrete frailty model
##'
##' This function plots estimates of the survival or cumulative hazard for each group, coloured according to the latent population that each group belongs to.
##'
##' @param x A fitted nonparametric discrete frailty model, as returned by \code{\link{npdf_cox}} with \code{estK=FALSE}.
##'
##' @param type character. If \code{"km"} group-specific Kaplan-Meier estimates of survival are plotted.   If \code{"na"} group-specific Nelson-Aalen estimates of the cumulative hazard are plotted.  The default is \code{"km"}.
##'
##' @param cols Vector of colour names or numbers, of the same length as the number of groups.  If not given, this defaults to \code{x$belonging}.
##'
##' @param survfit_opts Optional list of additional arguments to pass to \code{\link{survfit.formula}}.
##'
##' @param na_opts Optional list of arguments (other than \code{"cols"}) to pass to \code{\link{plot.nelsonaalen_npdf}}.
##'
##' @param ... Optional arguments to pass to \code{\link{plot.survfit}}.
##'
##' @examples
##'
##' result = npdf_cox( Surv(time, status) ~ x, groups = family, data = weibdata2030,
##'                     K = 2, estK = FALSE, eps_conv=10^-4)
##' plot( result )
##' plot( result, type = "km" )
##' plot( result, cols = ifelse( result$belonging == 1, "purple", "black" ), xlim = c( 0, 150 ) )
##'
##' ## use of survfit_opts.  show only first 10 groups
##' plot( result, survfit_opts = list(subset = (weibdata2030$family >= 10) ))
##'
##' plot( result, type = "na" )
##'
##' ## use of na_opts to customise the Nelson-Aalen plot
##' plot( result, type = "na", cols=ifelse(result$belonging==2, "purple", "black"),
##'      na_opts = list(xlim=c(0,200), ylim=c(0,2),
##'                     xlab="Follow-up days",
##'                     ylab="Nelson-Aalen estimate"))
##'
##' @importFrom graphics lines plot lines.default plot.default
##' @importFrom stats knots stepfun plot.stepfun
##'
##' @export
plot.npdf <- function(x, type="km", cols=NULL, survfit_opts = NULL, na_opts = NULL, ...){
    type <- match.arg(type, c("km", "na"))
    if (is.null(cols)){
      cols = rep( 0, length( x$belonging ) )
      for( i in unique( x$belonging ) )
      {
        cols[ x$belonging == i ] = max( x$belonging ) - i + 1
      }
    }
    if( type=="na"){
        object <- c(list(nelsonaalen_npdf(x)), list(cols=cols))
        do.call("plot.nelsonaalen_npdf", c(object, na_opts))
    } else{
        object <- survfit_npdf(x, survfit_opts)
        plot.survfit_npdf(object, cols=cols, ...)
    }
    invisible()
}

##' Survival or cumulative hazard curves from a nonparametric discrete frailty model chosen by automatic model selection
##'
##' @param x An object returned by \code{\link{npdf_cox}} with \code{estK=TRUE}, containing a list of fitted nonparametric discrete frailty models
##'
##' @param K The number of latent populations which identifies the specific model to plot estimates from.  By default this is the model selected by the criterion specified when calling \code{\link{npdf_cox}}.
##'
##' @param ... Options to customise the plot, passed to \code{\link{plot.npdf}}.  See \code{\link{plot.npdf}} for a description of these.
##'
##' @examples
##' result = npdf_cox( Surv(time, status) ~ x, groups = family, data = weibdata2030,
##'                    K = 2, eps_conv=10^-4)
##' plot( result, K = 2 )
##' plot( result, type = "na" )
##' plot( result, type = "na", cols=ifelse(result$belonging==2, "purple", "black"),
##'      na_opts = list(xlim=c(0,200), ylim=c(0,2),
##'                     xlab="Follow-up days",
##'                     ylab="Nelson-Aalen estimate"))
##'
##' @seealso \code{\link{plot.npdf}}
##'
##' @export
plot.npdflist <- function(x, K=NULL, ...){
    if (is.null(K)) K <- x$Kopt[x$criterion]
    mod <- x$models[[K]]
    mod$mf <- x$mf
    plot(mod, ...)
    invisible()
}

##' Nelson-Aalen estimates of group-specific cumulative hazard from a nonparametric discrete frailty model
##'
##' @inheritParams plot.npdf
##'
##' @return A list of objects of length equal to the number of groups in the data.  Each component is a list, equivalent to the output of \code{\link{survfit}} called for the corresponding group with \code{type="fh"}, but with two additional components:
##'
##' \code{y0}:  \code{-log} of the survival estimate
##'
##' \code{sfun0}:  a step function for \code{y0} in terms of time, estimated using \code{\link{stepfun}}.
##'
##' @seealso \code{\link{plot.nelsonaalen_npdf}}, \code{\link{plot.npdf}}.
##'
##' @examples
##' x = npdf_cox( Surv(time, status) ~ x, groups=family, data=weibdata2030, K = 2,
##'                estK=FALSE, eps_conv=10^-4)
##' object = nelsonaalen_npdf( x )
##'
##' @export
nelsonaalen_npdf <- function(x){
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
    class(res) <- "nelsonaalen_npdf"
    res
}

##' Plot Nelson-Aalen estimates of group-specific cumulative hazard from a nonparametric discrete frailty model
##'
##' @param x Object returned by \code{\link{nelsonaalen_npdf}} representing Nelson-Aalen estimates from a nonparametric discrete frailty model
##'
##' @param xlim x-axis limits (vector of 2)
##'
##' @param ylim x-axis limits (vector of 2)
##'
##' @param xlab x-axis label
##'
##' @param ylab y-axis label
##'
##' @param cols vector of colour names or numbers, of the same length as the number of groups
##'
##' @param ... options to pass to the generic \code{plot} function
##'
##' @examples
##'
##'  x = npdf_cox( Surv(time, status) ~ x, groups=family, data=weibdata2030, K = 2,
##'                 estK=FALSE, eps_conv=10^-4)
##' object = nelsonaalen_npdf( x )
##' plot( object )
##' plot( object, xlim=c(0,200), ylim=c(0,2),
##'      xlab="Follow-up days", ylab="Nelson-Aalen estimate",
##'      cols=ifelse(x$belonging==2, "purple", "black"))
##'
##' @seealso \code{\link{nelsonaalen_npdf}}
##'
##' @export
plot.nelsonaalen_npdf <- function(x, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, cols=NULL, ...){
    ngrp <- length(x)
    sfuns <- lapply(x, "[[", "sfun0")
    if (is.null(xlim)) xlim <- range(lapply(sfuns, knots))
    if (is.null(ylim)) ylim <- c(0,1)
    if (is.null(xlab)) xlab <- 'Time'
    if (is.null(ylab)) ylab <- 'Cumulative hazard'
    if (is.null(cols)) cols <- attr(x, "belonging")
    plot.default( NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
    for (i in seq(length=ngrp)){
        plot.stepfun( sfuns[[i]], col = cols[ i ], add=TRUE, do.points=FALSE )
    }
}

##' Kaplan-Meier estimates of group-specific cumulative hazard from a nonparametric discrete frailty model
##'
##' @inheritParams plot.npdf
##'
##' @return A list of survival estimates, one for each group, as produced by \code{\link{survfit.formula}}
##'
##' @seealso \code{\link{plot.npdf}}, \code{\link{plot.survfit_npdf}}
##'
##' @export
survfit_npdf <- function(x, survfit_opts = NULL){
    Y <- model.extract(x$mf, "response")
    Y <- as.data.frame(unclass(Y))
    groups <- model.extract(x$mf, "groups")
    survfit.call <- list(formula = Surv(time, status) ~ groups,
                         data = Y)
    res <- do.call("survfit", c(survfit.call, survfit_opts))
    class(res) <- "survfit_npdf"
    res
}

##' Plot Kaplan-Meier estimates of group-specific cumulative hazard from a nonparametric discrete frailty model
##'
##' @param x Object returned by \code{\link{survfit_npdf}} representing Kaplan-Meier estimates from a nonparametric discrete frailty model
##'
##' @inheritParams plot.nelsonaalen_npdf
##'
##' @seealso \code{\link{plot.npdf}}, \code{\link{survfit_npdf}}
##'
##' @examples
##'  x = npdf_cox( Surv(time, status) ~ x, groups=family, data=weibdata2030, K = 2,
##'                 estK=FALSE, eps_conv=10^-4)
##' object = survfit_npdf( x )
##' plot( object )
##'
##' @importFrom utils getFromNamespace
##'
##' @export
plot.survfit_npdf <- function(x, xlab=NULL, ylab=NULL, cols=NULL, ...){
    # plot.survfit chooses sensible defaults for axis limits
    if (is.null(xlab)) xlab <- 'Time'
    if (is.null(ylab)) ylab <- 'Survival'
    if (is.null(cols)) cols <- attr(x, "belonging")
    class(x) <- "survfit"
    plot(x, xlab=xlab, ylab=ylab, col=cols, ...)
}
