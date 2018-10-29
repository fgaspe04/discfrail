##' discfrail: Cox models for time-to-event data with nonparametric discrete group-specific frailties.
##'
##' discfrail: Functions for fitting Cox proportional hazards models for grouped time-to-event data, where the shared group-specific frailties have a discrete nonparametric distribution. Groups of individuals are clustered into a number of latent "populations" of groups, each with a common frailty that is assumed to act multiplicatively on their hazard.   Covariates can be included in the manner of the standard Cox model.   The baseline hazard is left unspecified as in the Cox model, and estimation is performed by maximum partial likelihood through an EM algorithm. 
##'
##' There are also functions for simulating from these models, and from similar models with a parametric baseline survival function.
##'
##' The methods are fully described and illustrated in Gasperoni et al. (2018). 
##'
##' \code{\link{npdf}} fits nonparametric discrete frailty models.  
##'
##' The number of latent populations of groups can either be fixed, or estimated through an automated model comparison and selection procedure. 
##'
##' \code{\link{plot.npdf}} plots fitted survival or cumulative hazard curves from the fitted models, by group, coloured by the latent population the group belongs to. 
##'
##' \code{\link{sim.npdf}} TODO Describe and rename
##'
##' 
##' @name discfrail-package
##' @aliases discfrail-package discfrail
##' @docType package
##'
##' @references
##' F. Gasperoni, F. Ieva, A.-M. Paganoni, C. Jackson, L. Sharples. (2018) 
##' Nonparametric frailty Cox models for hierarchical time-to-event data. 
##' Biostatistics (in press).
