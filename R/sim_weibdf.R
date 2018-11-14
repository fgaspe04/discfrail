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
#' @inherit sim_npdf return
#'
#' @details The "proportional hazards" parameterisation of the Weibull distribution is used, with survivor function \eqn{S(t) = exp(-lambda w exp(x^T beta) t^rho)}. Note this is different from the "accelerated failure time" parameterisation used in, e.g. \code{\link{dweibull}}.  Distribution functions for the proportional hazards parameterisation can be found in the \pkg{flexsurv} package.
#' 
#' @export
#'
#' @examples
#' N <- 100
#' S <- 40
#' lambda <- 0.5
#' beta <- 1.6
#' rho <- 1
#' p <- c( 0.8, 0.2 )
#' w_values <- c( 0.8, 1.6 )
#' data <- sim_weibdf( N, S, lambda, rho, beta, p, w_values)
#' head( data )
#'

sim_weibdf <- function( N, S = NULL, lambda, rho, beta, p, w_values )
{
  # if S is NULL, we sample the clusters' size from a Poisson with mean = 50
  if( is.null(S) ){
    S <- rpois( N, 50  )
  }

  # n is the total sample size
  n <- ifelse( length( S ) > 1, sum( S ), S*N )

  # covariate is sample from a normal
  x <- matrix( 0, nrow = n, ncol = length( beta ) )

  for( i in 1:length(beta))
  {
    x[ ,i] <- rnorm( n, 0, 1)
  }

  # frailty term for each group
  w <- rep( w_values, round( p*N ) )

  # because of round function we can have length( w ) != N
  if( length( w ) < N )
  {
    w <- c( w, tail( w, 1 ) )
  }else if( length( w ) > N )
  {
    length( w ) <-  N
  }

  #we shuffle the w
  index_shuffled = sample( 1:N )
  w = w[ index_shuffled ]

  # if S is given as a number, all groups are of the same size
  if( length( S ) == 1 ){
    coef = rep( S, N )
  }else{
    coef = S
  }

  # frailty term for each individual
  w_tot = rep( w, coef  )

  # computing event times
  v <- runif( n )
  Tlat <- (- log(v) / (lambda * w_tot * exp(x %*% beta) ) )^(1 / rho)

  # fixing parameters for the censoring distribution
  mean_cens <- mean(Tlat)
  var_cens <- mean(Tlat)/10

  # computing censoring times
  C <- rnorm( n, mean_cens, var_cens)

  # follow-up times and event indicators
  time <- pmin( Tlat, C )
  status <- as.numeric( Tlat <= C )

  # data set
  data.frame(family=rep( 1:N, coef ),
             time=time,
             status=status,
             x=x,
             belong = w_tot)
}
