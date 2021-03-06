N <- 100
S <- 40
Lambda_0_inv = function( t, c=0.01, d=4.6 ) ( t^( 1/d ) )/c
beta <- 1.6
p <- c( 0.8, 0.2 )
w_values <- c( 0.8, 1.6 )
set.seed(1)
weibdata <- sim_npdf( N, S, beta, Lambda_0_inv, p, w_values)

N <- 20 
S <- 30
Lambda_0_inv = function( t, c=0.01, d=4.6 ) ( t^( 1/d ) )/c
beta <- 1.6
p <- c( 0.6, 0.4 )
w_values <- c( 0.8, 3.2 )
set.seed(1)
weibdata20 <- sim_npdf( N, S, beta, Lambda_0_inv, p, w_values)

library(devtools)
use_data(weibdata)
use_data(weibdata20)
