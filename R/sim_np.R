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
#' @param N number of groups in the data
#' @param S number of individuals in each group
#' @param beta vector of log hazard ratios
#' @param Lambda_0_inv inverse cumulative baseline hazard function, that is, with covariate values 0 and frailty ratio 1
#' @param p vector of K elements.  The kth element gives the proportion of groups in the kth latent population of groups.
#' @param w_values vector of K distinct frailty values, one for each latent population.
#'
#' @return A data frame with one row for each simulated individual, and the following columns: 
#'
#' \code{family}: the group which the individual is in (integers 1, 2, ...)
#'
#' \code{time}: the simulated event time 
#'
#' \code{status}: the simulated survival status.   Censoring times are generated from a normal distribution with mean given by the simulated event time, and variance given by the mean divided by 10.  The event time is observed if it is less than the censoring time, and censored otherwise.
#'
#' \code{x}: matrix of covariate values, generated from a standard normal distribution.
#'
#' \code{belong}:  the frailty hazard ratio corresponding to the the cluster of groups in which the individual's group has been allocated. 
#'
#'
#' @export
#'
#' @examples
#' N <- 100
#' S <- 40
#' Lambda_0_inv = function( t, c=0.01, d=4.6 ) ( t^( 1/d ) )/c
#' beta <- 1.6
#' p <- c( 0.8, 0.2 )
#' w_values <- c( 0.8, 1.6 )
#' data <- simulWeibDiscreteFrailCovNPbaseInv( N, S, beta, Lambda_0_inv, p, w_values)
#' head( data )


simulWeibDiscreteFrailCovNPbaseInv <- function( N, S = NULL, beta, Lambda_0_inv, p, w_values )
{
  # if S is NULL, we sample the clusters' size from a Poisson with mean = 50
  if( is.null(S) )
    S <- rpois( N, 50  )

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
  if( length(w) < N )
  {
    w <- c( w, tail( w, 1 ) )
  }else if( length( w ) > N )
  {
    length( w ) <- N
  }

  #we shuffle the w
  index_shuffled <- sample( 1:N )
  w <- w[ index_shuffled ]

  # if S is given as a number, all groups are of the same size
  if( length( S ) == 1 ){
    coef <- rep( S, N )
  }else{
    coef <- S
  }

  # frailty term for each individual
  w_tot <- rep( w, coef )

  # # computing event times with Lambda_0_inv
  v <- runif( n )
  Tlat <- Lambda_0_inv(- log(v) / ( w_tot * exp(x %*% beta) ) )

  # fixing parameters for the censoring distribution
  mean_cens <- mean(Tlat)
  var_cens <- mean(Tlat)/10

  # censoring times
  C <- rnorm( n , mean_cens, var_cens )

  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)

  # data set
  data.frame(family=rep( 1:N, coef ),
             time=time,
             status=status,
             x=x,
             belong = w_tot)
}
