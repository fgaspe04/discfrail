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
        Y$pop <- rep(belonging, table(groups))
        print( Y )
        print( table( Y$pop ) )
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
