#' Generating time-to-event data with a nonparametric baseline
#'
#' @param N number of groups/clusters in data
#' @param S number of statistical units in each group
#' @param beta vector of regression parameters
#' @param Lambda_0_inv inverse cumulative baseline function
#' @param p vector of proportions related to latent populations
#' @param w_values vector of frailty values
#'
#' @return This fuction returns a dataset which fits a semiparametric Cox model with a multiplicative discrete frailty term.
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
