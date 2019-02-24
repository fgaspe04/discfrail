# Copyright (C) 2017 Francesca Gasperoni <francesca.gasperoni@polimi.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Simulation of grouped time-to-event data with Weibull baseline hazard and discrete shared frailty distribution
#'
#' This function returns a dataset generated from a Weibull proportional hazards model with a shared discrete frailty term, for given Weibull parameters, hazard ratios, distribution of groups among latent populations, frailty values for each latent population, and randomly-generated covariate values.
#'
#' @inheritParams sim_npdf
#'
#' @param lambda Weibull baseline rate parameter (see below), interpreted as the rate parameter with covariate values of 0 and frailty ratio 1.  For \eqn{rho=1} this is the baseline hazard.
#'
#' @param rho Weibull shape parameter (see below)
#'
#' @param beta covariate effects in the Weibull distribution, interpreted as log hazard ratios (see below)
#'
#' @param p vector of K elements.  The kth element gives the proportion of groups in the kth latent population of groups.
#'
#' @param w_values vector of K distinct frailty values, one for each latent population.
#'
#' @param cens_perc cens_perc percentage of censored events. Censoring times are assumed to be distributed as a Normal with variance equal to 1.
#'
#' @inherit sim_npdf return
#'
#' @details The "proportional hazards" parameterisation of the Weibull distribution is used, with survivor function \eqn{S(t) = exp(-\lambda t^{\rho} w exp(x^T {\beta}) )}. Note this is different from the "accelerated failure time" parameterisation used in, e.g. \code{\link{dweibull}}.  Distribution functions for the proportional hazards parameterisation can be found in the \pkg{flexsurv} package.
#'
#' @references
#' Wan, F. (2017). Simulating survival data with predefined censoring rates for proportional hazards models. \emph{Statistics in medicine}, 36(5), 838-854.
#'
#' @importFrom stats dnorm qnorm ecdf integrate uniroot
#'
#' @export
#'
#' @examples
#' J <- 100
#' N <- 40
#' lambda <- 0.5
#' beta <- 1.6
#' rho <- 1
#' p <- c( 0.8, 0.2 )
#' w_values <- c( 0.8, 1.6 )
#' cens_perc <- 0.2
#' data <- sim_weibdf( J, N, lambda, rho, beta, p, w_values, cens_perc)
#' head( data )
#'

sim_weibdf <- function( J, N = NULL, lambda, rho, beta, p, w_values, cens_perc )
{
  # if N is NULL, we sample the clusters' size from a Poisson with mean = 50
  if( is.null(N) ){
    N <- rpois( J, 50  )
  }

  # n is the total sample size
  n <- ifelse( length( N ) > 1, sum( N ), N*J )

  # covariate is sample from a normal
  x <- matrix( 0, nrow = n, ncol = length( beta ) )

  for( i in 1:length(beta))
  {
    x[ ,i] <- rnorm( n, 0, 1)
  }

  # frailty term for each group
  w <- rep( w_values, round( p*J ) )

  # because of round function we can have length( w ) != J
  if( length( w ) < J )
  {
    w <- c( w, tail( w, 1 ) )
  }else if( length( w ) > J )
  {
    length( w ) <-  J
  }

  #we shuffle the w
  index_shuffled = sample( 1:J )
  w = w[ index_shuffled ]

  # if N is given as a number, all groups are of the same size
  if( length( N ) == 1 ){
    coef = rep( N, J )
  }else{
    coef = N
  }

  # frailty term for each individual
  w_tot = rep( w, coef  )

  # computing event times
  v <- runif( n )
  Tlat <- (- log(v) / (lambda * w_tot * exp(x %*% beta) ) )^(1 / rho)

  # check cens_perc
  stopifnot(0 <= cens_perc && cens_perc <= 1)

  #
  esurv = ecdf(Tlat)
  efail = function(t) 1 - esurv(t)
  censor.mu = uniroot(function(mu)
    cens_perc -
      (integrate(function(t) efail(t)*dnorm(t, mu, 1),
                 -Inf, qnorm(0.5, mu, 1), subdivisions=2000 )$value +
         integrate(function(t) efail(t)*dnorm(t, mu, 1),
                   qnorm(0.5, mu, 1), Inf, subdivisions=2000 )$value),
    lower=0, upper=100, extendInt="up")$root

  Censor = rnorm( n, censor.mu, 1)

  # follow-up times and event indicators
  time <- pmin( Tlat, Censor )
  status <- as.numeric( Tlat <= Censor )

  # data set
  data.frame(family=rep( 1:J, coef ),
             time=time,
             status=status,
             x=x,
             belong = w_tot)
}
