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

#' Simulation of grouped time-to-event data with nonparametric baseline hazard and discrete shared frailty distribution
#'
#' This function returns a dataset generated from a semiparametric proportional hazards model with a shared discrete frailty term, for given cumulative baseline hazard function, hazard ratios, distribution of groups among latent populations, frailty values for each latent population, and randomly-generated covariate values.
#'
#' @param J number of groups in the data
#'
#' @param N number of individuals in each group
#'
#' @param beta vector of log hazard ratios
#'
#' @param Lambda_0_inv inverse cumulative baseline hazard function, that is, with covariate values 0 and frailty ratio 1
#'
#' @param p vector of K elements.  The kth element gives the proportion of groups in the kth latent population of groups.
#'
#' @param w_values vector of K distinct frailty values, one for each latent population.
#'
#' @param cens_perc percentage of censored events. Censoring times are assumed to be distributed as a Normal with variance equal to 1.
#'
#' @return A data frame with one row for each simulated individual, and the following columns:
#'
#' \code{family}: the group which the individual is in (integers 1, 2, ...)
#'
#' \code{time}: the simulated event time.
#'
#' \code{status}: the simulated survival status. Censoring times are generated from a Normal distribution with standard deviation equal to 1 and the mean is estimated in order to guarantee the determined percentage of censored events. The event time is observed (status=1) if it is less than the censoring time, and censored otherwise (status=0).
#'
#' \code{x}: matrix of covariate values, generated from a standard normal distribution.
#'
#' \code{belong}:  the frailty hazard ratio corresponding to the cluster of groups in which the individual's group has been allocated.
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
#' Lambda_0_inv = function( t, c=0.01, d=4.6 ) ( t^( 1/d ) )/c
#' beta <- 1.6
#' p <- c( 0.8, 0.2 )
#' w_values <- c( 0.8, 1.6 )
#' cens_perc <- 0.1
#' data <- sim_npdf( J, N, beta, Lambda_0_inv, p, w_values, cens_perc )
#' head( data )


sim_npdf <- function( J, N = NULL, beta, Lambda_0_inv, p, w_values, cens_perc )
{
  # if N is NULL, we sample the clusters' size from a Poisson with mean = 50
  if( is.null(N) )
    N <- rpois( J, 50  )

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
  if( length(w) < J )
  {
    w <- c( w, tail( w, 1 ) )
  }else if( length( w ) > J )
  {
    length( w ) <- J
  }

  #we shuffle the w
  index_shuffled <- sample( 1:J )
  w <- w[ index_shuffled ]

  # if N is given as a number, all groups are of the same size
  if( length( N ) == 1 ){
    coef <- rep( N, J )
  }else{
    coef <- N
  }

  # frailty term for each individual
  w_tot <- rep( w, coef )

  # computing event times with Lambda_0_inv
  v <- runif( n )
  Tlat <- Lambda_0_inv(- log(v) / ( w_tot * exp(x %*% beta) ) )


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
  time <- pmin(Tlat, Censor)
  status <- as.numeric(Tlat <= Censor)

  # data set
  data.frame(family=rep( 1:J, coef ),
             time=time,
             status=status,
             x=x,
             belong = w_tot)
}
