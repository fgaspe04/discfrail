#' Hierarchical time-to-event data simulated from a Weibull baseline distribution with a nonparametric frailty
#' 
#' @format A data frame with 4000 rows and the following columns 
#'
#' \describe{
#' \item{family}{Integer from 1 to 100 indicating the group membership}
#' \item{time}{Survival time}
#' \item{status}{Indicator for an observed event at the survival time}
#' \item{x}{Covariate value}
#' \item{belong}{Underlying group-specific frailty value}
#' }
#'
#' TODO document values used for simulation 
#' 
"weibdata"

#' Hierarchical time-to-event data simulated from a Weibull baseline distribution with a nonparametric frailty
#' 
#' @format A data frame with 600 rows and the following columns 
#'
#' \describe{
#' \item{family}{Integer from 1 to 20 indicating the group membership}
#' \item{time}{Survival time}
#' \item{status}{Indicator for an observed event at the survival time}
#' \item{x}{Covariate value}
#' \item{belong}{Underlying group-specific frailty value}
#' }
#'
#' TODO document values used for simulation 
#' 
"weibdata20"
