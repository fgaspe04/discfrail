###################################################
#######################################################
#
#     Summary of all functions updated
#
#######################################################
#######################################################

rm(list=ls())
library(survival)
library(numDeriv)
library(Matrix)
#library(xtable)

simulWeibDiscreteFrail = function(N, S = NULL, lambda, rho, beta, p, w_values)
{
  #set.seed(1234)
  if( is.null(S) )
    S = rpois( N, 50  )

  #n is the total sample size
  n = ifelse( length( S ) > 1, sum(S), S )

  # covariate --> random from a normal
  x = rnorm( n, 0, 1)

  # Frailty term
  w = rep( w_values, round( p*N ) )

  if( length(w) < N )
  {
    w = c(w,tail(w,1))
  }
  else if(length(w) > N)
  {
    length( w ) = N
  }

  #set.seed(1234)
  index_shuffled = sample( 1:N )

  #try to shuffle the w
  w = w[ index_shuffled ]

  coef = S

  if( length( S ) == 1 )
    coef = rep( S, N )

  w_tot = rep( w, coef  )

  # Weibull latent event times
  v <- runif( n )
  Tlat <- (- log(v) / (lambda * w_tot * exp( x * beta ) ) )^(1 / rho)

  mean_cens = mean(Tlat)
  var_cens = mean(Tlat)/10

  # censoring times
  C <- rnorm( n, mean_cens, var_cens)

  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)

  # data set
  data.frame(family=rep( 1:N, coef ),
             #rep = rep( 1:S, N ),
             time=time,
             status=status,
             x=x,
             belong = w_tot)
}

simulWeibDiscreteFrailCovNPbaseInv <- function(N, S=NULL, beta, Lambda_0_inv, p, w_values)
{
  #set.seed(1234)
  if( is.null(S) )
    S = rpois( N, 50  )

  #n is the total sample size
  n = ifelse( length( S ) > 1, sum(S), S*N )

  x = matrix(0, nrow=n, ncol=length(beta))
  # covariate --> random from a normal
  for( i in 1:length(beta))
  {
    x[ ,i] = rnorm( n, 0, 1)
  }

  # Frailty term
  w = rep( w_values, round( p*N ) )

  if( length(w) < N )
  {
    w = c(w,tail(w,1))
  }
  if(length(w) > N)
  {
    length( w ) = N
  }

  #set.seed(1234)
  index_shuffled = sample( 1:N )

  #try to shuffle the w
  w = w[ index_shuffled ]

  coef = S

  if( length( S ) == 1 )
    coef = rep( S, N )

  w_tot = rep( w, coef )

  # Lambda_0_inv
  v <- runif( n )
  Tlat <- Lambda_0_inv(- log(v) / ( w_tot * exp(x %*% beta) ) )
  #Lambda_0 <- Vectorize(function(y) uniroot(function(x)
  #  Lambda_0_inv(x) - y, lower=0, upper=100, extendInt="yes")$root)

  mean_cens = mean(Tlat)
  var_cens = mean(Tlat)/10

  # censoring times
  C <- rnorm( n , mean_cens, var_cens)


  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)

  # data set
  data.frame(family=rep( 1:N, coef ),
             #rep = rep(1:S, N ),
             time=time,
             status=status,
             x=x,
             belong = w_tot)
}


##########################
# Computing the s.e matrix
##########################

exact_sec_der = function( ncovs, w, p, D, sumhaz, dsumhaz, d2sumhaz, lamsum)
{
  K = length(w)
  H = length(D)

  #to build
  d2p = array(dim=c(H, K-1, K-1))
  d2w = array(dim=c(H, K, K))
  d2cov = array(dim=c(H, ncovs, ncovs))

  d2pw = array(dim=c(H, K-1, K))
  d2pcov = array(dim=c(H, K-1, ncovs))
  d2wcov = array(dim=c(H, K, ncovs))

  exact_der = matrix(0, nrow=(K-1+K+ncovs), ncol=(K-1+K+ncovs))


  if( K > 1 )
  {
    w_sum_haz = matrix( unlist( lapply( 1:K, function(i) w[i] * sumhaz )), nrow = H, ncol= K)
    denominator = rowSums(matrix( unlist( lapply( 1:K, function(i) p[i]*w[i]^D )), nrow = H, ncol= K)*exp(-w_sum_haz))
    numerator1 = rowSums(matrix( unlist( lapply( 1:K, function(i) p[i]*w[i]^(D+1) )), nrow = H, ncol= K)*exp(-w_sum_haz))
    numerator2 = rowSums(matrix( unlist( lapply( 1:K, function(i) p[i]*w[i]^(D+2) )), nrow = H, ncol= K)*exp(-w_sum_haz))
  }else{
    denominator = numerator1 = numerator2 = exp(-sumhaz)
  }
  check=0

  for( j in 1:H)
  {

    if(K>=2){
      for( g in 1:(K-1) )
      {
        #Pure term 1 #dim(dpdp) = (K-1)x(K-1)
        for( l in 1:(K-1) )
        {
          d2p[ j, g, l] = -(w[g]^D[j]*exp(-w[g]*sumhaz[j])-w[K]^D[j]*exp(-w[K]*sumhaz[j]))*(w[l]^D[j]*exp(-w[l]*sumhaz[j])-w[K]^D[j]*exp(-w[K]*sumhaz[j]))/(denominator[j])^2

        }

        #Mixed term 1-2 #dim(d2pw) = (K-1)xK
        for(l in 1:K)
        {
          d2pw[ j, g, l ] = -(w[g]^D[j]*exp(-w[g]*sumhaz[j])-w[K]^D[j]*exp(-w[K]*sumhaz[j]))*(p[l]*w[l]^(D[j]-1)*exp(-w[l]*sumhaz[j])*(D[j]-w[l]*sumhaz[j]))/(denominator[j])^2+
            + ifelse( l == g, w[l]^D[j]*exp(-w[l]*sumhaz[j])*(D[j]-w[l]*sumhaz[j])/denominator[j], 0) +
            + ifelse( l == K, w[l]^D[j]*exp(-w[l]*sumhaz[j])*(D[j]-w[l]*sumhaz[j])/denominator[j], 0)
        }

        #Mixed term 1-3 #dim(d2pcov) = (K-1)xncovs
        for(r in 1:ncovs)
        {
          d2pcov[ j, g, r ] = dsumhaz[j,r] * ( (w[K]^(D[j]+1)*exp(-w[K]*sumhaz[j]) - w[g]^(D[j]+1)*exp(-w[g]*sumhaz[j]))/denominator[j] +
                                                 -(w[K]^(D[j])*exp(-w[K]*sumhaz[j]) - w[g]^(D[j])*exp(-w[g]*sumhaz[j]))*numerator1[j]/(denominator[j])^2 )
        }

      }
    }

    #Pure term 2 #dim(d2w) = KxK
    for( g in 1:K )
    {
      for( l in 1:K )
      {

        d2w[ j, g, l] = -p[g]*p[l]*((w[g]*w[l])^(D[j]-1))*exp(-(w[g]+w[l])*sumhaz[j])*(D[j]-w[g]*sumhaz[j])*(D[j]-w[l]*sumhaz[j])/(denominator[j])^2+
          + ifelse(g!=l, 0, p[g]*w[g]^(D[j]-2)*exp(-w[g]*sumhaz[j])*((D[j]-1-w[g]*sumhaz[j])*(D[j]-w[g]*sumhaz[j])-w[g]*sumhaz[j])/denominator[j])

      }
    }

    #Pure term 3 #dim(d2cov) = ncovs x ncovs
    for(r in 1:ncovs)
    {
      for(s in 1:ncovs)
      {
        d2cov[ j, r, s ] = lamsum[ j, r, s ]+ numerator2[j]/denominator[j]*(dsumhaz[j,r]*dsumhaz[j,s])+
          - numerator1[j]/denominator[j]*d2sumhaz[j,r,s] +
          - (numerator1[j]/denominator[j])^2*(dsumhaz[j,r]*dsumhaz[j,s])
      }
    }

    #Mixed term 2-3 #dim(d2wcov) = Kxncovs
    for(g in 1:K)
    {
      for(r in 1:ncovs)
      {
        d2wcov[ j, g, r ] = p[g] * w[g]^(D[j]-1) * exp(-w[g]*sumhaz[j]) * dsumhaz[j,r] * ( (-w[g]*(D[j]+1-w[g]*sumhaz[j]))/denominator[j]+
                                                                                             + ((D[j]-w[g]*sumhaz[j])*numerator1[j])/(denominator[j])^2)
      }
    }

    check = ifelse( length(which(is.na(d2w[j,,])))>0,j,check)
    if(K>2){
      exact_dertmp = rbind( cbind(d2p[j, , ], d2pw[j, , ], d2pcov[j, , ]),
                            cbind(t(d2pw[j, , ]), d2w[j, , ], d2wcov[j, , ]),
                            cbind(t(d2pcov[j, , ]), t(d2wcov[j, , ]), d2cov[j, , ]) )
    }else if( K==2 ){
      exact_dertmp = rbind( c(d2p[j, , ], d2pw[j, , ], d2pcov[j, , ]),
                            cbind(d2pw[j, , ], d2w[j, , ], d2wcov[j, , ]),
                            cbind(d2pcov[j, , ], t(d2wcov[j, , ]), d2cov[j, , ]) )
    }else if(K==1){
      exact_dertmp = rbind( c(d2w[j, , ], d2wcov[j, , ]), cbind( d2wcov[j, , ], d2cov[j, , ] ) )
    }

    exact_der = exact_der  + exact_dertmp

  }
  return(exact_der)

}

exact_first_der = function(ncovs, w, p, D, sumhaz, dsumhaz, lamroutlam0, deltaijcovr, dsumhazred )
{
  K = length(w)
  H = length(D)

  #to build
  dp = rep(0,K-1)
  dw = rep(0,K)
  dcov = rep(0,ncovs)
  dcov2 = rep(0,ncovs)

  exact_derivative1 = exact_derivative2 = rep(0, K-1+K+ncovs)

  if( K > 1 )
  {
    w_sum_haz = matrix( unlist( lapply( 1:K, function(i) w[i] * sumhaz )), nrow = H, ncol= K)
    denominator = rowSums(matrix( unlist( lapply( 1:K, function(i) p[i]*w[i]^D )), nrow = H, ncol= K)*exp(-w_sum_haz))
    numerator1 = rowSums(matrix( unlist( lapply( 1:K, function(i) p[i]*w[i]^(D+1) )), nrow = H, ncol= K)*exp(-w_sum_haz))
    numerator2 = rowSums(matrix( unlist( lapply( 1:K, function(i) p[i]*w[i]^(D+2) )), nrow = H, ncol= K)*exp(-w_sum_haz))
  }else{
    denominator = numerator1 = numerator2 = exp(-sumhaz)
  }


  for( j in 1:H )
  {
    for(g in 1:(K-1))
    {
      dp[g] = (w[g]^(D[j])*exp(-w[g]*sumhaz[j]) - w[K]^(D[j])*exp(-w[K]*sumhaz[j]))/denominator[j]
    }

    for(g in 1:K)
    {
      dw[g] = (p[g] * w[g]^(D[j]-1) * exp(-w[g]*sumhaz[j]) * (D[j]-w[g]*sumhaz[j]))/denominator[j]
    }

    for(r in 1:ncovs)
    {
      dcov[r] = -numerator1[j]*dsumhaz[j,r]/denominator[j] + lamroutlam0[j,r]
      dcov2[r] = -numerator1[j]*dsumhazred[j,r]/denominator[j] + deltaijcovr[j,r]
    }

    exact_derivative1 = exact_derivative1 + c( dp, dw, dcov)
    exact_derivative2 = exact_derivative2 + c( dp, dw, dcov2)
  }

  return(list(exact_derivative1,exact_derivative2))

}




Iy = function( ncovs, alpha, w, p, D, sumhaz, lamroutlam0, dsumhaz)#, lamwoutlam0, dsumhazw)
{
  #ncovs: number of covariates
  #alpha: matrix of posterior density of Z, HxK
  #w: vector of frailties, Kx1
  #p: vector of p, Kx1
  #D: number of events in each group Hx1
  #sumhaz: Lambda_0*exp(X_ij* beta) H x 1
  #lamroutlam0: lambda_0r/lambda_0 + X_ijr  H x nocovs
  #dsumhaz: first derivative of cumulative baseline haz H x ncovs

  K = dim(alpha)[2]
  H = dim(alpha)[1]

  #to build
  dpdp = array(dim=c(H, K-1))
  dwdw = array(dim=c(H, K))
  dcovdcov = array(dim=c(H, ncovs))

  expectedStexpectedSj = matrix(0, nrow=(K-1+K+ncovs), ncol=(K-1+K+ncovs))

  #each array is group specific

  for( j in 1:H)
  {
    #Pure term 1 #dim(dpdp) = (K-1)xq
    dpdp[ j, ] =  alpha[ j, 1:(K-1) ]/p[1:(K-1)] - alpha[ j, K ]/p[K]

    #Pure term 2 #dim(dwdw) = K x K

    dwdw[ j, ] = alpha[ j, 1:K]*(D[j]/w[1:K]-sumhaz[j]) #( lamwoutlam0[ j, 1:K ] - dsumhazw[j,1:K] )


    #Pure term 3 #dim(dcovdcov) = ncovs x ncovs
    for(r in 1:ncovs) {
      dcovdcov[ j, r ] = alpha[ j, ] %*% t(lamroutlam0[j,r]-w*dsumhaz[j,r])
    }


    expectedStexpectedSj = expectedStexpectedSj + c( dpdp[j,], dwdw[j, ], dcovdcov[j, ] ) %*% t(c( dpdp[j,], dwdw[j, ], dcovdcov[j, ] ))


  }

  return(expectedStexpectedSj)
}


Ixy = function( ncovs, alpha, w, p, D, sumhaz, lamroutlam0, dsumhaz ) #, lamwoutlam0, dsumhazw)
{
  #ncovs: number of covariates
  #alpha: matrix of posterior density of Z, HxK
  #w: vector of frailties, Kx1
  #p: vector of p, Kx1
  #D: number of events in each group Hx1
  #sumhaz: Lambda_0*exp(X_ij* beta) H x 1
  #lamroutlam0: lambda_0r/lambda_0 + X_ijr  H x nocovs
  #dsumhaz: first derivative of cumulative baseline haz H x ncovs

  K = dim(alpha)[2]
  H = dim(alpha)[1]

  #to build
  dpdp = array(dim=c(H, K-1, K-1))
  dwdw = array(dim=c(H, K, K))
  dcovdcov = array(dim=c(H, ncovs, ncovs))

  dpdw = array(dim=c(H, K-1, K))
  dpdcov = array(dim=c(H, K-1, ncovs))
  dwdcov = array(dim=c(H, K, ncovs))

  expectedStSgroupj = matrix(0, nrow=(K-1+K+ncovs), ncol=(K-1+K+ncovs))

  #each array is group specific

  for( j in 1:H)
  {
    #Pure term 1 #dim(dpdp) = (K-1)x(K-1)
    dpdp[ j, , ] = alpha[ j, K ]/p[K]^2
    if( K > 2){
      diag( dpdp[ j, , ] ) = diag( dpdp[ j, , ] ) + alpha[ j, 1:(K-1) ]/p[1:(K-1)]^2
    }else if( K <= 2){
      dpdp[ j, , ] = alpha[ j, K ]/p[K]^2 + alpha[ j, 1:(K-1) ]/p[1:(K-1)]^2
    }


    #Pure term 2 #dim(dwdw) = K x K
    #diag( alpha[ j, 1:K]*( D[j]/w[1:K] - sumhaz[j] )^2 )
    if(K>1){
      dwdw[ j, , ] = diag(alpha[ j, 1:K]*(D[j]/w[1:K]-sumhaz[j])^2)#( lamwoutlam0[ j, 1:K ] - dsumhazw[j,1:K] )^2 )
    }else if(K==1){
      dwdw[ j, , ] = alpha[ j, 1:K]*(D[j]/w[1:K]-sumhaz[j])^2#( lamwoutlam0[ j, 1:K ] - dsumhazw[j,1:K] )^2 )
    }

    #Pure term 3 #dim(dcovdcov) = ncovs x ncovs
    for(r in 1:ncovs) {
      for(s in 1:ncovs) {
        dcovdcov[ j, r, s ] = alpha[j,] %*% t( (lamroutlam0[j,r]-w*dsumhaz[j,r])*(lamroutlam0[j,s]-w*dsumhaz[j,s]) )
      }
    }

    #Mixed term 1-2 #dim(dpdw) = (K-1)xK
    dpdwtmp = rep(0,(K-1))
    dpdwtmp[K-1] = -alpha[ j,K]/p[K] * (D[j]/w[K]-sumhaz[j])
    if(K>2){
      dpdw[j, , ] = cbind( diag(alpha[ j,1:(K-1)]/p[1:(K-1)] * (D[j]/w[1:(K-1)]-sumhaz[j])), dpdwtmp )
    }else{
      dpdw[j, , ] = cbind( alpha[ j,1:(K-1)]/p[1:(K-1)] * (D[j]/w[1:(K-1)]-sumhaz[j]), dpdwtmp )
    }

    #Mixed term 1-3 #dim(dpdcov) = K-1 x ncovs
    for (r in 1:ncovs)
    {
      dpdcov[j, , r] = alpha[j,1:(K-1)]/p[1:(K-1)] * (lamroutlam0[j,r]-w[1:(K-1)]*dsumhaz[j,r]) - alpha[j,K]/p[K]*(lamroutlam0[j,r]-w[K]*dsumhaz[j,r])
    }



    #Mixed term 2-3 #dim(dpdcov) = K x ncovs
    for (r in 1:ncovs)
    {
      dwdcov[j, , r] = alpha[j,1:K] * (lamroutlam0[j,r]-w[1:K]*dsumhaz[j,r]) * (D[j]/w[1:K]-sumhaz[j])
    }

    if(K>2){
      expectedStSgroupjtmp = rbind( cbind(dpdp[j, , ], dpdw[j, , ], dpdcov[j, , ]),
                                    cbind(t(dpdw[j, , ]), dwdw[j, , ], dwdcov[j, , ]),
                                    cbind(t(dpdcov[j, , ]), t(dwdcov[j, , ]), dcovdcov[j, , ]) )
    }else if( K==2 ){
      expectedStSgroupjtmp = rbind( c(dpdp[j, , ], dpdw[j, , ], dpdcov[j, , ]),
                                    cbind(dpdw[j, , ], dwdw[j, , ], dwdcov[j, , ]),
                                    cbind(dpdcov[j, , ], t(dwdcov[j, , ]), dcovdcov[j, , ]) )
    }else if(K==1){
      expectedStSgroupjtmp = rbind( c(dwdw[j, , ], dwdcov[j, , ]), cbind( dwdcov[j, , ], dcovdcov[j, , ] ) )
    }

    expectedStSgroupj = expectedStSgroupj  + expectedStSgroupjtmp

  }

  return(expectedStSgroupj)
}

Ix = function( ncovs, alpha, w, p, D, d2sumhaz, lamsum, dsumhaz) #,lamw2outlam02, d2sumhazw, lamwbetaoutlam02 ,d2sumhazwbeta )
{
  #ncovs: number of covariates
  #alpha: matrix of posterior density of Z, HxK
  #w: vector of frailties, Kx1
  #p: vector of p, Kx1
  #D: number of events in each group Hx1
  #d2sumhaz: second derivative of cumulative baseline haz H x ncovs x ncovs
  #lamsum: second derivative of baseline haz H x nocovs x ncovs
  #dsumhaz: first derivative of cumulative baseline haz H x ncovs

  K = dim(alpha)[2]
  H = dim(alpha)[1]

  #pure terms
  d2pp = array( dim = c( H, K-1, K-1 ) )
  d2ww = array( dim = c( H, K, K ) )
  d2covcov = array( dim = c( H, ncovs, ncovs ) )

  #mixed terms
  d2wcov = array( dim = c( H, K , ncovs ) )

  #result
  expectedBgroupj = matrix(0, nrow = (K-1+K+ncovs), ncol = (K-1+K+ncovs) )

  for(j in 1:H)
  {
    #Pure term 1 #dim(d2pp) = (K-1)x(K-1)
    d2pp[ j, , ] = - alpha[ j, K ]/p[K]^2
    if(K>2){
      diag( d2pp[ j, , ] ) = diag( d2pp[ j, , ] ) - alpha[ j, 1:(K-1) ]/p[1:(K-1)]^2
    }else{
      d2pp[ j, , ] = - alpha[ j, K ]/p[K]^2 - alpha[ j, 1:(K-1) ]/p[1:(K-1)]^2
    }


    #Pure term 2 #dim(d2pp) = K x K
    if(K>1){
      d2ww[ j, , ] = diag( -alpha[ j, 1:K ]*D[j]/w[1:K]^2 ) #alpha[ j, 1:K ]*(lamw2outlam02[j,1:K]-d2sumhazw[j,1:K] )
    }else if(K==1){
      d2ww[ j, , ] = -alpha[ j, 1:K ]*D[j]/w[1:K]^2 #alpha[ j, 1:K ]*(lamw2outlam02[j,1:K]-d2sumhazw[j,1:K] )
    }

    #Pure term 3 #dim(d2covcov) = ncovs x ncovs
    for(r in 1:ncovs)
    {
      for(s in 1:ncovs)
      {
        d2covcov[ j, r, s ] = t(alpha[ j, 1:K]) %*% (lamsum[ j, r, s]- w[1:K] * d2sumhaz[ j , r, s])
      }
    }

    #Mixed term 1 #dim(d2wcov) = K x ncovs
    for(k in 1:K)
    {
      for(r in 1:ncovs)
      {
        d2wcov[ j, k, r ] = - alpha[j, k] * dsumhaz[ j, r ] #(lamwbetaoutlam02[j, k, r]-d2sumhazwbeta[j, k, r])
      }
    }

    if(K>=2)
    {
      expectedBgroupjtmp = rbind( cbind( d2pp[ j, , ], matrix( 0, nrow = K-1, ncol = K ), matrix( 0, nrow = K-1, ncol = ncovs ) ),
                                  cbind( matrix( 0, nrow = K, ncol = K-1 ), d2ww[j, , ], d2wcov[j, , ] ),
                                  cbind( matrix( 0, nrow = ncovs, ncol = K-1 ), t(d2wcov[j, , ]), d2covcov[j, , ]) )
    }else{
      expectedBgroupjtmp = rbind( c( d2ww[j, , ], d2wcov[j, , ] ),
                                  cbind( d2wcov[j, , ], d2covcov[j, , ]) )
    }
    expectedBgroupj = expectedBgroupj + expectedBgroupjtmp

  }

  return(expectedBgroupj)
}

information_matrix = function( expectedS, expectedB, expectedSexpectedS )
{
  info = -expectedB-expectedS+expectedSexpectedS
  cov <- solve(info)
  se <- sqrt(diag(cov))

  return(list(cov,se))
}



npdf_cox <- function( data, time, status, family, covariates, population = 2, eps_conv=10^-4)
{
  H <- length( unique( data[[family]] ) ) #total number of hospital
  K <- population
  dat_ord <- data[ order( data[[family]], data[[time]], -data[[status]] ), ]
  nt <- length(unique(data[[time]]))
  ncovs <- length(covariates)
  cov_read <- paste(covariates, sep=',')
  cum_haz <- as.data.frame( cbind( hazard=rep( 0, nt ),
                                   time = sort( unique( data[[time]] ) )) )
  sumhaz <- rep( 0, H )
  dsumhaz <- array( dim = c( H, ncovs ) ) #perchè non ho fatto una matrix??
  d2sumhaz <- lamsum <- array( dim = c( H, ncovs, ncovs ) )
  deltaijcovr <- matrix( 0, nrow = H, ncol = ncovs )
  dsumhazred <- matrix( 0, nrow = H, ncol = ncovs )

  #### computing SE
  #cum_haz$hazard <-
  YY <- haz <- rep( 0, nt )
  dYY <- dhaz <- dcumhaz <- matrix( nrow = nt, ncol = ncovs )
  #dYYmatr <- dhazmatr <- matrix( nrow = nt, ncol = ncovs )

  d2YY <- d2haz <- d2cumhaz <- array( dim = c( nt, ncovs, ncovs ) )
  #d2YYmatr <- d2hazmatr <- array( dim = c( nt, ncovs, ncovs ) )

  # #new line
  lamroutlam0 <- matrix( nrow = H, ncol = ncovs )

  ##sort the cluster id (there is a smarter way)
  for( i in 1:length(unique( dat_ord[[family]] ) ) )
  {
    value <- unique( dat_ord[[family]] )[i]
    dat_ord[[family]][ dat_ord[[family]] == value ] <- i
  }
  #
  D <- table( factor( dat_ord[[family]] )[ dat_ord[[status]] == 1 ] )  # number of events in each group j
  n_H <- as.vector( table( dat_ord[[family]] ) ) #number of patients in each hospital

  risk_index <- matrix( 0, nrow = dim(dat_ord)[1], ncol = nt ) # matrix of at risk patients: nrows x nt
  N = rep( 0, nt )

  time_list <- sapply( 1:nt, function(x) !is.na(match(dat_ord[[time]], cum_haz$time[x]))) #matrix dim(dat_ord)[1] x nt, time_list[i,j] = T if the i-th row in the data is the j-th time recorded
  N <- sapply( 1:dim(time_list)[2], function(x) sum(dat_ord[[status]][time_list[,x]])) #number of events happenend at each time. N[j] is the number of eents happenend at time j


  for( l in 1:nt )
  {
    risk_index[ which( dat_ord[[time]] >= cum_haz$time[ l ]), l ] <-  1 #nrow(dat_ord) x nt
  }

  lis = list()

  while( K >= 1 )
  {
    print(K)
    count <- 0
    not_proper_K <- 0
    eps <- 10^5

    p <- rep( 1/K, K )
    w <- seq( 1, 5, length.out = K )

    if( K==1 ) w <- 1

    set.seed(2)
    w_to_rep <- sample( rep( w, round( p*H ) ) )
    if( H - sum( round( p*H ) ) >= 1 )
    {
      w_to_rep <- c( w_to_rep, rep( w_to_rep[ sum( round( p*H ) ) ], H - sum( round( p*H ) ) ) )
      p <- table(w_to_rep)/H
    }
    if( H - sum( round( p*H ) ) <= -1 )
    {
      length( w_to_rep ) <- H
      p <- table(w_to_rep)/H
    }

    w_off <- w_to_rep[ dat_ord[[family]] ]

    numerator <- rep( 0, K )
    E_formula_part <- rep( 0, H )
    alpha <- E_formula <- matrix( 0, nrow = H, ncol = K)


    while(eps > eps_conv & count < 200 )
    {

      mF <- formula( paste("Surv(time, status)", paste( paste( covariates, collapse= '+'), "offset( log( w_off ) )", sep ='+'), sep = "~"))
      temp_model <- survival::coxph( mF, data = dat_ord,  method = "breslow" )


      beta_hat <- temp_model$coef

      # computing baseline hazard and cumulative baseline hazard
      YY = as.numeric( (w_off*exp( beta_hat %*% t(dat_ord[,cov_read]) )) %*% risk_index ) # nt x 1 contribute from at risk patient
      haz = N/YY
      cum_haz$hazard = cumsum( haz )

      # saving p
      p_old <- p


      #Expectation step

      for( j in 1:H )
      {
        hospj <- dat_ord[[family]]==j
        ebz <- exp( as.matrix( dat_ord[hospj,cov_read] ) %*% beta_hat )
        ti <- match(dat_ord[[time]][hospj], cum_haz$time)
        lam0t <- cum_haz$hazard[ti]


        # flag_cov is an indicator for the presence of at least 1 covariate
        flag_cov = ifelse( length(cov_read) != 0, 1, 0 )

        E_formula_part[j] <- ifelse( flag_cov == 1,
                                     sum( lam0t*ebz ),
                                     sum( lam0t ) )

        for( k in 1:K )
        {
          E_formula[j,k] <- ifelse( flag_cov == 1,
                                    sum( lam0t*ebz*w[ k ] ),
                                    sum( lam0t*w[ k ] ) )


          numerator[ k ] <- p[ k ] * exp( D[ j ] * log( w[ k ] ) +
                                            -E_formula[j,k] )

        }

        alpha[ j, ] <- numerator/sum( numerator )

      }

      #Maximization Step

      # vector of latent population label, dim = H x 1
      belonging <- as.numeric( apply( alpha, 1, which.max) )


      # vector of proportions
      p <- ( colSums( alpha ) )/ H
      print(p)

      # vector of frailty, dim = Kx1
      w <- (D %*% alpha)/(E_formula_part %*% alpha)

      # standardization of w
      w <- w/min(w)
      #w <- w/w[1]
      #w[which.min(w)] <- 1

      # vector of frailty values, dim = nx1
      w_off <- ( w[ belonging ][dat_ord[[family]]] )

      # computation of the distance from the previous estimate of p
      eps <- max( abs(p-p_old) )

      # augmenting count
      count <- count + 1
      print(count)

    }

    # computing 1st and 2nd hazard derivatives
    dhazards <- harzard_derivatives( dat_ord, risk_index, w_off, beta_hat, cov_read, N, YY )
    dhaz <- dhazards$dhaz
    d2haz <- dhazards$d2haz
    dcumhaz <- dhazards$dcumhaz
    d2cumhaz <- dhazards$d2cumhaz

    print('dhaz done')

    #nt x ncovs
    # dYY <- t( risk_index) %*% ( as.matrix(sapply(w_off * exp( beta_hat %*% t(dat_ord[ ,cov_read ]) ) * dat_ord[ ,cov_read ], as.numeric)))
    # dhaz <- - dYY * N/YY^2
    #
    # if( ncovs > 1 )
    # {
    #   #nrow(dat_ord) x ncovs x ncovs
    #   d2YYtemp <- array( t( apply( dat_ord[,cov_read], 1, function(x) x%*%t(x)) ), dim=c( dim( dat_ord )[1], ncovs, ncovs ) )
    #   #nt x ncovs x ncovs
    #   dYYdYY <- array( t( apply( dYY[,1:ncovs], 1, function(x) x%*%t(x)) ), dim=c( nt, ncovs, ncovs ) )
    #   L <- lapply( seq_len( dim(risk_index)[2] ), function(i) apply(as.vector(w_off*exp(beta_hat %*% t(dat_ord[,cov_read])))*d2YYtemp*risk_index[,i],MARGIN = 2:3,sum))
    #   d2YY <- array( unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)))
    #   d2YY <- aperm( d2YY, c( 3, 1, 2 ))
    # }else{
    #   d2YYtemp <- dat_ord[,cov_read]^2
    #   dYYdYY <- dYY^2
    #   d2YY <- t(risk_index) %*% ( w_off*exp( beta_hat * dat_ord[,cov_read] )*d2YYtemp ) # nt x ncovs
    # }
    #
    # # constructing array for point-wise operation
    # d2haz <- array( N/YY^3*(2*dYYdYY-d2YY*YY), dim = c( nt, ncovs, ncovs ))
    #
    #
    # if(ncovs>1)
    # {
    #   dcumhaz <- matrix( apply( dhaz, MARGIN = 2, cumsum ), ncol = ncovs) #sum by column
    #   d2cumhaz <- array( apply( d2haz, MARGIN = 2:3, cumsum ), dim = c( nt, ncovs, ncovs ) )
    # }else{
    #   dcumhaz <- matrix( cumsum( dhaz ), ncol = ncovs)
    #   d2cumhaz <- array( cumsum( d2haz ), dim = c( nt, ncovs, ncovs ))
    # }

    #computing the input for the se functions
    input_se <- input_se_functions( family, time, status, dat_ord, cov_read, beta_hat, haz, cum_haz, dhaz, dcumhaz, d2haz, d2cumhaz )
    lamroutlam0 <- input_se$lamroutlam0
    lamsum <- input_se$lamsum
    deltaijcovr <- input_se$deltaijcovr
    sumhaz <- input_se$sumhaz
    dsumhazred <- input_se$dsumhazred
    dsumhaz <- input_se$dsumhaz
    d2sumhaz <- input_se$d2sumhaz

    print('imput_se done')

    # for( j in 1:H )
    # {
    #   hospj <- dat_ord[[family]]==j
    #   ebz <- exp( as.matrix( dat_ord[hospj,cov_read] ) %*% beta_hat )
    #   ti <- match(dat_ord[[time]][hospj], cum_haz$time)
    #   lam0t <- cum_haz$hazard[ti]
    #   sumhaz[j] <- sum(lam0t * ebz)
    #
    #
    #   for (r in 1:ncovs){
    #     cov_int_r = cov_read[r]
    #     covr <- dat_ord[hospj,cov_int_r]
    #     lamrt <- dcumhaz[ti, r]
    #     #dsumhaz is the second part of III term
    #     dsumhaz[j,r] <- sum(ebz*(covr*lam0t + lamrt))
    #     #new lines for III term
    #     lamroutlam0tmp <- dat_ord[[status]][hospj]*( dhaz[ti,r]/haz[ti] + covr )
    #     lamroutlam0[j,r] <- sum(lamroutlam0tmp[haz[ti]>0])
    #
    #     deltaijcovrtmp <- dat_ord[[status]][hospj]*covr
    #     deltaijcovr[j,r] <- sum(deltaijcovrtmp[haz[ti]>0])
    #     dsumhazred[j,r] <- sum(ebz*(covr*lam0t))
    #
    #     for (s in 1:ncovs){
    #       cov_int_s <- cov_read[s]
    #       covs <- dat_ord[hospj,cov_int_s]
    #       lamst <- dcumhaz[ti, s]
    #       lamrst <- d2cumhaz[ti, r, s]
    #       d2sumhaz[j,r,s] <- sum(ebz*(covr*covs*lam0t + covr*lamst + covs*lamrt + lamrst))
    #       lamtmp <- dat_ord[[status]][hospj] * (d2haz[ti,r,s]*haz[ti] - dhaz[ti,r]*dhaz[ti,s])/haz[ti]^2
    #       lamsum[j,r,s] <- sum(lamtmp[haz[ti]>0]) # exclude zero hazard
    #
    #     }
    #
    #   }
    # }

    #function to be numerically optimized
    sefundeponwbeta = function( x )
    {
      YYfun <- as.numeric( (w_off*exp( x[(2*K):(2*K-1+ncovs)] %*% t(dat_ord[,cov_read]) )) %*% risk_index ) # nt x 1 contribute from at risk patient
      hazfun <- N/YYfun
      cum_hazfun <- cumsum(hazfun)

      lfull <- rep(0,H)
      for(j in 1:H)
      {
        hospj <- dat_ord[[family]]==j
        ebz <- exp( as.matrix( dat_ord[hospj,cov_read] ) %*% x[(2*K):(2*K-1+ncovs)] )
        ti <- match(dat_ord[[time]][hospj], cum_haz$time)
        lam0t <- cum_hazfun[ti]
        haz0t <- hazfun[ti]


        firsttermtmp <- dat_ord[[status]][hospj]*log(haz0t*ebz)
        firstterm <- sum(firsttermtmp[haz0t>0,])

        secondtermtmp <- -sum(lam0t * ebz) * x[(K):(2*K-1)]
        wpowerD <- x[(K):(2*K-1)]^D[j]
        secondterm <- log(sum(exp(secondtermtmp)*wpowerD*c(x[1:(K-1)],1-sum(x[1:(K-1)]))))

        lfull[j] <- firstterm+secondterm
      }
      return( sum( lfull ) )
    }

    #Numerical computation of standard errors
    if( length( unique( belonging ) ) == K )
    {
      if( K != 1 )
      {
        deriv <- numDeriv::genD( sefundeponwbeta, c( p[1:(K-1)], w, beta_hat ), method = "Richardson" )
      }else{
        deriv <- numDeriv::genD( sefundeponwbeta, c( w, beta_hat), method = "Richardson" )
      }
      info <- matrix(0, 2*K-1+ncovs, 2*K-1+ncovs)
      info[upper.tri(info, diag=T)] <- deriv$D[(2*K+ncovs):(2*K-1+ncovs+(2*K+ncovs-1)*(2*K+ncovs)/2)]
      info <- Matrix::forceSymmetric(info)
      infoRich <- solve(-info)

      print('infoRich done')

      #Louis method for standard errors
      EStS <- Ixy( ncovs, alpha, w, p, D, sumhaz, lamroutlam0, dsumhaz)
      EB <- Ix( ncovs, alpha, w, p, D, d2sumhaz, lamsum, dsumhaz)
      ESES <- Iy( ncovs, alpha, w, p, D, sumhaz, lamroutlam0, dsumhaz)
      stderr_covariance <- information_matrix( expectedS = EStS, expectedB = EB, expectedSexpectedS = ESES )
      print('Louis done')

      #Exact second derivativefor standard errors
      SecDer <- exact_sec_der( ncovs, w, p, D, sumhaz, dsumhaz, d2sumhaz,lamsum )
      InfoFischer <- solve( -SecDer )
      print('InfoFischer done')
    }


    # computing loglikelihood, AIC and BIC
    llik <- log_likelihood( family, time, status, dat_ord, cov_read, p, w, alpha, beta_hat, haz, cum_haz)

    # llik = 0
    #
    # for( j in 1:H )
    # {
    #   hospj <- dat_ord[[family]]==j
    #   hospj_lam <- dat_ord[[family]]==j & dat_ord[[status]]==1
    #   ti <- match(dat_ord[[time]][hospj], cum_haz$time)
    #   ti_lam <- match(dat_ord[[time]][hospj_lam], cum_haz$time)
    #
    #   for( k in 1:K )
    #   {
    #     llik_LAM <- ifelse( flag_cov == 1,
    #                         sum( cum_haz$hazard[ti]*w[k]*exp( as.matrix( dat_ord[hospj,cov_read] ) %*% beta_hat ) ),
    #                         sum( cum_haz$hazard[ti]*w[k] ) )
    #
    #     #w_off[hospj][1] because w_off[hospj] are all equal
    #     llik_lam <- ifelse( flag_cov == 1,
    #                         sum( log( haz[ti_lam]*w[k]*exp( as.matrix( dat_ord[hospj_lam,cov_read] ) %*% beta_hat ) ) ),
    #                         sum( log( haz[ti_lam]*w[k] ) ) )
    #
    #     llik_p <- ifelse( K==1, 1, p[ k ] ) #/( 1-p[ belonging[ j ] ] ) )
    #
    #     llik <- llik + alpha[j,k]*(llik_lam-llik_LAM)  + alpha[j,k]*log( llik_p )#*length(hospj)
    #
    #   }
    # }


    tot_param <- K-1 + K + length(cov_read)
    BIC <- -2*llik + tot_param * log( length( which(dat_ord$status == 1 ) ) )
    AIC <- -2*llik + 2*tot_param


    if( length( unique( belonging ) ) < K )
    {
      not_proper_K <- 1
      stderr_covariance <- NULL
      infoRich <- NULL
      InfoFischer <- NULL
    }

    step = K

    lis[[step]] <- list()
    lis[[step]][['K']] <- K
    lis[[step]][['too_much_K']] <- not_proper_K
    lis[[step]][['p']] <- p
    lis[[step]][['w']] <- w
    lis[[step]][['beta']] <- beta_hat
    lis[[step]][['belonging']] <- belonging
    lis[[step]][['llik']] <- llik
    lis[[step]][['BIC']] <- BIC
    lis[[step]][['AIC']] <- AIC
    lis[[step]][['varcovLouis']] <- stderr_covariance[[1]]
    lis[[step]][['seLouis']] <- stderr_covariance[[2]]
    lis[[step]][['varcovRich']] <- infoRich
    lis[[step]][['seRich']] <- sqrt(diag(infoRich))
    lis[[step]][['varcovExact']] <- InfoFischer
    lis[[step]][['seExact']] <- sqrt(diag(InfoFischer))


    K = K-1
  }

  best_llik = lis[[1]]$llik
  best_BIC = lis[[1]]$BIC
  best_AIC = lis[[1]]$AIC
  pos_llik = 1
  pos_BIC = 1
  pos_AIC = 1
  max_K = 0
  right_K = 1

  #Finding the maximum K for which we have distinct w
  for( i in 2:population)
  {
    if( lis[[i]]$too_much_K > max_K )
    {
      max_K = lis[[i]]$too_much_K
      right_K = lis[[i-1]]$K
    }
  }

  if(max_K==0)
  {
    right_K = lis[[population]]$K
  }
  #I use the right K as the top value for evaluating the best BIC and AIC
  i=1
  while( i <= right_K)
  {
    if( lis[[i]]$llik > best_llik)
    {
      best_llik = lis[[i]]$llik
      pos_llik = lis[[i]]$K
    }
    if( lis[[i]]$BIC < best_BIC)
    {
      best_BIC = lis[[i]]$BIC
      pos_BIC = lis[[i]]$K
    }
    if( lis[[i]]$AIC < best_AIC)
    {
      best_AIC = lis[[i]]$AIC
      pos_AIC = lis[[i]]$K
    }
    i=i+1
  }
  return( c(lis,pos_llik,pos_BIC,pos_AIC, right_K) ) #pos is the right K up to BIC or AIC, while right_K is the right K up to the belonging
}


#
# ##test
# N=100
# #S=40
# #S=rep(c(40,60,40,60,50),20)
# lambda = 0.5
# beta = c(1.6,3) #2#c(1.5,2,3)#c(0.05,0.26,0.36)
# rho=1
# p = c(0.8,0.2)#c(0.25,0.35,0.3,0.1)
# w_values= c(0.8,1.6)#c(1,2,2.7,3.2)
#
#
# data = simulWeibDiscreteFrail( N, S=NULL, lambda, rho, beta, p, w_values)
# test_res = npdf_cox( data, time = 'time', status = 'status', family = 'family', covariates = c('x'), population = 3, eps_conv=10^-4)
#
# set.seed(2000)
# data = simulWeibDiscreteFrailCovNPbaseInv(N, S=50, beta, Lambda_0_inv = function(t, c=0.01, d=4.6) (t^(1/d))/c , p, w_values)
# start.time <- Sys.time()
# test_res = npdf_cox( data, time = 'time', status = 'status', family = 'family', covariates = c('x.1','x.2'), population = 3, eps_conv=10^-4) #.1','x.2','x.3
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
#
#
# dat_hand = sdatright_trans_3
# test_real_data = npdf_cox( dat_hand, time = 'time', status = 'status', family = 'OSP_ID', covariates = c('AGE','SEX'), population = 5, eps_conv=10^-4)
#
#
# test_res[[1]]$w
# test_res[[2]]$w
# test_res[[3]]$w
#
# test_res[[1]]$p
# test_res[[2]]$p
# test_res[[3]]$p
#
# test_res[[1]]$beta
# test_res[[2]]$beta
# test_res[[3]]$beta
#
# test_res[[1]]$seLouis
# test_res[[1]]$seRich
# test_res[[1]]$seExact
# test_res[[2]]$seLouis
# test_res[[2]]$seRich
# test_res[[2]]$seExact
# test_res[[3]]$seLouis
# test_res[[3]]$seRich
# test_res[[3]]$seExact
#
# test_res[[5]]
# test_res[[6]]
# test_res[[7]]
#
# plot(1:3,c(test_res[[1]]$AIC,test_res[[2]]$AIC,test_res[[3]]$AIC))
# table( rep( test_res[[2]]$belonging, S)- data$belong )
# (350+320)/sum(table( rep( test_res[[2]]$belonging, S)- data$belong ))
#
