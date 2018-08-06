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

# Computing standard errors
#
# @param ncovs number of covariates
# @param w: vector of frailties, Kx1
# @param p: vector of proportions, Kx1
# @param D: vector of death per group, Hx1
# @param sumhaz: cumulative baseline
# @param dsumhaz: first derivative of cumulative master
# @param d2sumhaz: second derivative of cumulative master
# @param lamsum:

exact_sec_der <- function( ncovs, w, p, D, sumhaz, dsumhaz, d2sumhaz, lamsum)
{
  K <- length( w )
  H <- length( D )

  #to build
  d2p <- array(dim=c(H, K-1, K-1))
  d2w <- array(dim=c(H, K, K))
  d2cov <- array(dim=c(H, ncovs, ncovs))

  d2pw <- array(dim=c(H, K-1, K))
  d2pcov <- array(dim=c(H, K-1, ncovs))
  d2wcov <- array(dim=c(H, K, ncovs))

  exact_der <- matrix(0, nrow=(K-1+K+ncovs), ncol=(K-1+K+ncovs))

  if( K > 1 )
  {
    w_sum_haz = matrix( unlist( lapply( 1:K, function(i) w[i] * sumhaz )), nrow = H, ncol= K)
    denominator = rowSums(exp( matrix( unlist( lapply( 1:K, function(i) ( log(p[i]) + D*log(w[i])) )), nrow = H, ncol= K) -w_sum_haz))
    numerator1 = rowSums(exp( matrix( unlist( lapply( 1:K, function(i) ( log(p[i]) + (D+1)*log(w[i])) )), nrow = H, ncol= K) -w_sum_haz))
    numerator2 = rowSums(exp( matrix( unlist( lapply( 1:K, function(i) ( log(p[i]) + (D+2)*log(w[i])) )), nrow = H, ncol= K) -w_sum_haz))
  }else{
    denominator = numerator1 = numerator2 = exp(-sumhaz)
  }

  for( j in 1:H)
  {

    if(K>=2){
      for( g in 1:(K-1) )
      {
        #Pure term 1 #dim(dpdp) = (K-1)x(K-1)
        for( l in 1:(K-1) )
        {
          d2p[ j, g, l] = -(exp(D[j]*log(w[g])-w[g]*sumhaz[j])-exp(D[j]*log(w[K])-w[K]*sumhaz[j]))*(exp(D[j]*log(w[l])-w[l]*sumhaz[j])-exp(D[j]*log(w[K])-w[K]*sumhaz[j]))/(denominator[j])^2
        }

        #Mixed term 1-2 #dim(d2pw) = (K-1)xK
        for(l in 1:K)
        {
          d2pw[ j, g, l ] = -(exp(D[j]*log(w[g])-w[g]*sumhaz[j])-exp(D[j]*log(w[K])-w[K]*sumhaz[j]))*(p[l]*exp((D[j]-1)*log(w[l])-w[l]*sumhaz[j])*(D[j]-w[l]*sumhaz[j]))/(denominator[j])^2+
            + ifelse( l == g, exp(D[j]*log(w[l])-w[l]*sumhaz[j])*(D[j]-w[l]*sumhaz[j])/denominator[j], 0) +
            + ifelse( l == K, exp(D[j]*log(w[l])-w[l]*sumhaz[j])*(D[j]-w[l]*sumhaz[j])/denominator[j], 0)
        }

        #Mixed term 1-3 #dim(d2pcov) = (K-1)xncovs
        for(r in 1:ncovs)
        {
          d2pcov[ j, g, r ] = dsumhaz[j,r] * ( (exp((D[j]+1)*log(w[K])-w[K]*sumhaz[j]) - exp((D[j]+1)*log(w[g])-w[g]*sumhaz[j]))/denominator[j] +
                                                 -(exp(D[j]*log(w[K])-w[K]*sumhaz[j]) - exp(D[j]*log(w[g])-w[g]*sumhaz[j]))*numerator1[j]/(denominator[j])^2 )
        }

      }
    }

    #Pure term 2 #dim(d2w) = KxK
    for( g in 1:K )
    {
      for( l in 1:K )
      {

        d2w[ j, g, l] = -p[g]*p[l]*exp((D[j]-1)*log(w[g]*w[l])-(w[g]+w[l])*sumhaz[j])*(D[j]-w[g]*sumhaz[j])*(D[j]-w[l]*sumhaz[j])/(denominator[j])^2+
          + ifelse(g!=l, 0, p[g]*exp((D[j]-2)*log(w[g])-w[g]*sumhaz[j])*((D[j]-1-w[g]*sumhaz[j])*(D[j]-w[g]*sumhaz[j])-w[g]*sumhaz[j])/denominator[j])

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
        d2wcov[ j, g, r ] = p[g] * exp( (D[j]-1)*log(w[g]) - w[g]*sumhaz[j]) * dsumhaz[j,r] * ( (-w[g]*(D[j]+1-w[g]*sumhaz[j]))/denominator[j]+
                                                                                                  + ((D[j]-w[g]*sumhaz[j])*numerator1[j])/(denominator[j])^2)
      }
    }


    if( K > 2 ){
      exact_dertmp = rbind( cbind(d2p[j, , ], d2pw[j, , ], d2pcov[j, , ]),
                            cbind(t(d2pw[j, , ]), d2w[j, , ], d2wcov[j, , ]),
                            cbind(t(d2pcov[j, , ]), t(d2wcov[j, , ]), d2cov[j, , ]) )
    }else if( K == 2 ){
      exact_dertmp = rbind( c(d2p[j, , ], d2pw[j, , ], d2pcov[j, , ]),
                            cbind(d2pw[j, , ], d2w[j, , ], d2wcov[j, , ]),
                            cbind(d2pcov[j, , ], t(d2wcov[j, , ]), d2cov[j, , ]) )
    }else if( K == 1 ){
      exact_dertmp = rbind( c(d2w[j, , ], d2wcov[j, , ]), cbind( d2wcov[j, , ], d2cov[j, , ] ) )
    }

    exact_der = exact_der  + exact_dertmp

  }
  return(exact_der)
}

# @param ncovs number of covariates
# @param w: vector of frailties, Kx1
# @param p: vector of proportions, Kx1
# @param D: vector of death per group, Hx1
# @param sumhaz: Lambda_0*exp(X_ij* beta), H x 1
# @param dsumhaz: first derivative of cumulative baseline haz, H x ncovs
# @param d2sumhaz: second derivative of cumulative baseline haz, H x ncovs
# @param lamroutlam0: lambda_0r/lambda_0 + X_ijr, H x ncovs
# @param deltaijcovr:
# @param dsumhazred:

exact_first_der <- function( ncovs, w, p, D, sumhaz, dsumhaz, lamroutlam0, deltaijcovr, dsumhazred )
{
  K <- length( w ) #number of latent populations
  H <- length( D ) #number of groups

  #to build
  dp <- rep( 0, K-1 )
  dw <- rep( 0, K )
  dcov <- rep( 0, ncovs )
  dcov2 <- rep( 0, ncovs )

  exact_derivative1 <- exact_derivative2 <- rep( 0, K-1+K+ncovs )

  if( K > 1 )
  {
    w_sum_haz <- matrix( unlist( lapply( 1:K, function(i) w[i] * sumhaz )), nrow = H, ncol= K)
    denominator <- rowSums(matrix( unlist( lapply( 1:K, function(i) p[i]*w[i]^D )), nrow = H, ncol= K)*exp(-w_sum_haz))
    numerator1 <- rowSums(matrix( unlist( lapply( 1:K, function(i) p[i]*w[i]^(D+1) )), nrow = H, ncol= K)*exp(-w_sum_haz))
    numerator2 <- rowSums(matrix( unlist( lapply( 1:K, function(i) p[i]*w[i]^(D+2) )), nrow = H, ncol= K)*exp(-w_sum_haz))
  }else{
    denominator <- numerator1 <- numerator2 <- exp(-sumhaz)
  }


  for( j in 1:H )
  {
    for(g in 1:(K-1))
    {
      dp[g] <- (w[g]^(D[j])*exp(-w[g]*sumhaz[j]) - w[K]^(D[j])*exp(-w[K]*sumhaz[j]))/denominator[j]
    }

    for(g in 1:K)
    {
      dw[g] <- (p[g] * w[g]^(D[j]-1) * exp(-w[g]*sumhaz[j]) * (D[j]-w[g]*sumhaz[j]))/denominator[j]
    }

    for(r in 1:ncovs)
    {
      dcov[r] <- -numerator1[j]*dsumhaz[j,r]/denominator[j] + lamroutlam0[j,r]
      dcov2[r] <- -numerator1[j]*dsumhazred[j,r]/denominator[j] + deltaijcovr[j,r]
    }

    exact_derivative1 <- exact_derivative1 + c( dp, dw, dcov)
    exact_derivative2 <- exact_derivative2 + c( dp, dw, dcov2)
  }

  return( list( exact_derivative1, exact_derivative2 ) )
}

# @param ncovs: number of covariates
# @param alpha: matrix of posterior density of Z, HxK
# @param w: vector of frailties, Kx1
# @param p: vector of p, Kx1
# @param D: number of events in each group, Hx1
# @param sumhaz: Lambda_0*exp(X_ij* beta), H x 1
# @param lamroutlam0: lambda_0r/lambda_0 + X_ijr,  H x ncovs
# @param dsumhaz: first derivative of cumulative baseline haz, H x ncovs

Iy <- function( ncovs, alpha, w, p, D, sumhaz, lamroutlam0, dsumhaz)
{
  K <- dim(alpha)[2] #number of latent populations
  H <- dim(alpha)[1] #number of groups

  #to build
  dpdp <- array(dim=c(H, K-1))
  dwdw <- array(dim=c(H, K))
  dcovdcov <- array(dim=c(H, ncovs))

  expectedStexpectedSj <- matrix(0, nrow=(K-1+K+ncovs), ncol=(K-1+K+ncovs))

  #each array is group specific

  for( j in 1:H)
  {
    #Pure term 1 #dim(dpdp) = (K-1)xq
    dpdp[ j, ] <- alpha[ j, 1:(K-1) ]/p[1:(K-1)] - alpha[ j, K ]/p[K]

    #Pure term 2 #dim(dwdw) = K x K

    dwdw[ j, ] <- alpha[ j, 1:K]*(D[j]/w[1:K]-sumhaz[j])


    #Pure term 3 #dim(dcovdcov) = ncovs x ncovs
    for(r in 1:ncovs) {
      dcovdcov[ j, r ] <- alpha[ j, ] %*% t(lamroutlam0[j,r]-w*dsumhaz[j,r])
    }


    expectedStexpectedSj <- expectedStexpectedSj + c( dpdp[j,], dwdw[j, ], dcovdcov[j, ] ) %*% t(c( dpdp[j,], dwdw[j, ], dcovdcov[j, ] ))


  }

  return(expectedStexpectedSj)
}

# @param ncovs: number of covariates
# @param alpha: matrix of posterior density of Z, HxK
# @param w: vector of frailties, Kx1
# @param p: vector of p, Kx1
# @param D: number of events in each group, Hx1
# @param sumhaz: Lambda_0*exp(X_ij* beta), H x 1
# @param lamroutlam0: lambda_0r/lambda_0 + X_ijr,  H x ncovs
# @param dsumhaz: first derivative of cumulative baseline haz, H x ncovs

Ixy = function( ncovs, alpha, w, p, D, sumhaz, lamroutlam0, dsumhaz )
{
  K <- dim(alpha)[2] #number of latent populations
  H <- dim(alpha)[1] #number of groups

  #to build
  dpdp <- array(dim=c(H, K-1, K-1))
  dwdw <- array(dim=c(H, K, K))
  dcovdcov <- array(dim=c(H, ncovs, ncovs))

  dpdw <- array(dim=c(H, K-1, K))
  dpdcov <- array(dim=c(H, K-1, ncovs))
  dwdcov <- array(dim=c(H, K, ncovs))

  expectedStSgroupj <- matrix(0, nrow=(K-1+K+ncovs), ncol=(K-1+K+ncovs))

  for( j in 1:H)
  {
    #Pure term 1 #dim(dpdp) = (K-1)x(K-1)
    dpdp[ j, , ] <- alpha[ j, K ]/p[K]^2
    if( K > 2){
      diag( dpdp[ j, , ] ) <- diag( dpdp[ j, , ] ) + alpha[ j, 1:(K-1) ]/p[1:(K-1)]^2
    }else if( K <= 2){
      dpdp[ j, , ] <- alpha[ j, K ]/p[K]^2 + alpha[ j, 1:(K-1) ]/p[1:(K-1)]^2
    }


    #Pure term 2 #dim(dwdw) = K x K
    if(K>1){
      dwdw[ j, , ] <- diag(alpha[ j, 1:K]*(D[j]/w[1:K]-sumhaz[j])^2)#( lamwoutlam0[ j, 1:K ] - dsumhazw[j,1:K] )^2 )
    }else if(K==1){
      dwdw[ j, , ] <- alpha[ j, 1:K]*(D[j]/w[1:K]-sumhaz[j])^2#( lamwoutlam0[ j, 1:K ] - dsumhazw[j,1:K] )^2 )
    }

    #Pure term 3 #dim(dcovdcov) = ncovs x ncovs
    for(r in 1:ncovs) {
      for(s in 1:ncovs) {
        dcovdcov[ j, r, s ] <- alpha[j,] %*% t( (lamroutlam0[j,r]-w*dsumhaz[j,r])*(lamroutlam0[j,s]-w*dsumhaz[j,s]) )
      }
    }

    #Mixed term 1-2 #dim(dpdw) = (K-1)xK
    dpdwtmp <- rep(0,(K-1))
    dpdwtmp[K-1] <- -alpha[ j,K]/p[K] * (D[j]/w[K]-sumhaz[j])
    if(K>2){
      dpdw[j, , ] <- cbind( diag(alpha[ j,1:(K-1)]/p[1:(K-1)] * (D[j]/w[1:(K-1)]-sumhaz[j])), dpdwtmp )
    }else{
      dpdw[j, , ] <- cbind( alpha[ j,1:(K-1)]/p[1:(K-1)] * (D[j]/w[1:(K-1)]-sumhaz[j]), dpdwtmp )
    }

    #Mixed term 1-3 #dim(dpdcov) = K-1 x ncovs
    for (r in 1:ncovs)
    {
      dpdcov[j, , r] <- alpha[j,1:(K-1)]/p[1:(K-1)] * (lamroutlam0[j,r]-w[1:(K-1)]*dsumhaz[j,r]) - alpha[j,K]/p[K]*(lamroutlam0[j,r]-w[K]*dsumhaz[j,r])
    }

    #Mixed term 2-3 #dim(dpdcov) = K x ncovs
    for (r in 1:ncovs)
    {
      dwdcov[j, , r] <- alpha[j,1:K] * (lamroutlam0[j,r]-w[1:K]*dsumhaz[j,r]) * (D[j]/w[1:K]-sumhaz[j])
    }

    if(K>2){
      expectedStSgroupjtmp <- rbind( cbind(dpdp[j, , ], dpdw[j, , ], dpdcov[j, , ]),
                                    cbind(t(dpdw[j, , ]), dwdw[j, , ], dwdcov[j, , ]),
                                    cbind(t(dpdcov[j, , ]), t(dwdcov[j, , ]), dcovdcov[j, , ]) )
    }else if( K==2 ){
      expectedStSgroupjtmp <- rbind( c(dpdp[j, , ], dpdw[j, , ], dpdcov[j, , ]),
                                    cbind(dpdw[j, , ], dwdw[j, , ], dwdcov[j, , ]),
                                    cbind(dpdcov[j, , ], t(dwdcov[j, , ]), dcovdcov[j, , ]) )
    }else if(K==1){
      expectedStSgroupjtmp <- rbind( c(dwdw[j, , ], dwdcov[j, , ]), cbind( dwdcov[j, , ], dcovdcov[j, , ] ) )
    }

    expectedStSgroupj <- expectedStSgroupj  + expectedStSgroupjtmp

  }

  return(expectedStSgroupj)
}


# @param ncovs: number of covariates
# @param alpha: matrix of posterior density of Z, HxK
# @param w: vector of frailties, Kx1
# @param p: vector of proportion, Kx1
# @param D: number of events in each group, Hx1
# @param d2sumhaz: second derivative of cumulative baseline haz, H x ncovs x ncovs
# @param lamsum: second derivative of baseline haz, H x ncovs x ncovs
# @param dsumhaz: first derivative of cumulative baseline haz, H x ncovs

Ix <- function( ncovs, alpha, w, p, D, d2sumhaz, lamsum, dsumhaz)
{
  K <- dim(alpha)[2] #number of latent populations
  H <- dim(alpha)[1] #number of groups

  #pure terms
  d2pp <- array( dim = c( H, K-1, K-1 ) )
  d2ww <- array( dim = c( H, K, K ) )
  d2covcov <- array( dim = c( H, ncovs, ncovs ) )

  #mixed terms
  d2wcov <- array( dim = c( H, K , ncovs ) )

  #result
  expectedBgroupj <- matrix(0, nrow = (K-1+K+ncovs), ncol = (K-1+K+ncovs) )

  for(j in 1:H)
  {
    #Pure term 1 #dim(d2pp) = (K-1)x(K-1)
    d2pp[ j, , ] <- - alpha[ j, K ]/p[K]^2
    if(K>2){
      diag( d2pp[ j, , ] ) <- diag( d2pp[ j, , ] ) - alpha[ j, 1:(K-1) ]/p[1:(K-1)]^2
    }else{
      d2pp[ j, , ] <- - alpha[ j, K ]/p[K]^2 - alpha[ j, 1:(K-1) ]/p[1:(K-1)]^2
    }


    #Pure term 2 #dim(d2pp) = K x K
    if(K>1){
      d2ww[ j, , ] <- diag( -alpha[ j, 1:K ]*D[j]/w[1:K]^2 )
    }else if(K==1){
      d2ww[ j, , ] <- -alpha[ j, 1:K ]*D[j]/w[1:K]^2
    }

    #Pure term 3 #dim(d2covcov) = ncovs x ncovs
    for(r in 1:ncovs)
    {
      for(s in 1:ncovs)
      {
        d2covcov[ j, r, s ] <- t(alpha[ j, 1:K]) %*% (lamsum[ j, r, s]- w[1:K] * d2sumhaz[ j , r, s])
      }
    }

    #Mixed term 1 #dim(d2wcov) = K x ncovs
    for(k in 1:K)
    {
      for(r in 1:ncovs)
      {
        d2wcov[ j, k, r ] <- - alpha[j, k] * dsumhaz[ j, r ] #(lamwbetaoutlam02[j, k, r]-d2sumhazwbeta[j, k, r])
      }
    }

    if(K>=2)
    {
      expectedBgroupjtmp <- rbind( cbind( d2pp[ j, , ], matrix( 0, nrow = K-1, ncol = K ), matrix( 0, nrow = K-1, ncol = ncovs ) ),
                                  cbind( matrix( 0, nrow = K, ncol = K-1 ), d2ww[j, , ], d2wcov[j, , ] ),
                                  cbind( matrix( 0, nrow = ncovs, ncol = K-1 ), t(d2wcov[j, , ]), d2covcov[j, , ]) )
    }else{
      expectedBgroupjtmp <- rbind( c( d2ww[j, , ], d2wcov[j, , ] ),
                                  cbind( d2wcov[j, , ], d2covcov[j, , ]) )
    }
    expectedBgroupj <- expectedBgroupj + expectedBgroupjtmp

  }

  return( expectedBgroupj )
}

#computation of the information matrix

information_matrix <- function( expectedS, expectedB, expectedSexpectedS )
{
  info <- -expectedB-expectedS+expectedSexpectedS
  cov <- solve(info)
  se <- sqrt(diag(cov))

  return( list( cov, se ) )
}


hazard_derivatives <- function( X, risk_index, w_off, beta_hat, N, YY)
{
  ncovs <- ncol(X)
  nt <- dim(risk_index)[2]

  dhaz <- dcumhaz <- matrix( nrow = nt, ncol = ncovs )
  d2haz <- d2cumhaz <- array( dim = c( nt, ncovs, ncovs ) )

  #nt x ncovs
  hazbase <- as.numeric(exp( beta_hat %*% t(X) ))
  dYY <- t(risk_index) %*% (w_off * hazbase * X)
  dhaz <- - dYY * N/YY^2

  if( ncovs > 1 )
  {
      nobs <- nrow(X)
    d2YYtemp <- array( t( apply( X, 1, function(x) x%*%t(x)) ), dim=c( nobs, ncovs, ncovs ) )
    #nt x ncovs x ncovs
    dYYdYY <- array( t( apply( dYY[,1:ncovs], 1, function(x) x%*%t(x)) ), dim=c( nt, ncovs, ncovs ) )
    L <- lapply( seq_len( dim(risk_index)[2] ), function(i) apply(as.vector(w_off*exp(beta_hat %*% t(X)))*d2YYtemp*risk_index[,i],MARGIN = 2:3,sum))
    d2YY <- array( unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)))
    d2YY <- aperm( d2YY, c( 3, 1, 2 ))
  }else{
    d2YYtemp <- X^2
    dYYdYY <- dYY^2
    d2YY <- t(risk_index) %*% ( w_off*exp( beta_hat * X )*d2YYtemp ) # nt x ncovs
  }

  # constructing array for point-wise operation
  d2haz <- array( N/YY^3*(2*dYYdYY-d2YY*YY), dim = c( nt, ncovs, ncovs ))


  if(ncovs>1)
  {
    dcumhaz <- matrix( apply( dhaz, MARGIN = 2, cumsum ), ncol = ncovs) #sum by column
    d2cumhaz <- array( apply( d2haz, MARGIN = 2:3, cumsum ), dim = c( nt, ncovs, ncovs ) )
  }else{
    dcumhaz <- matrix( cumsum( dhaz ), ncol = ncovs)
    d2cumhaz <- array( cumsum( d2haz ), dim = c( nt, ncovs, ncovs ))
  }
  lis <- list()
  lis[['dhaz']] <- dhaz
  lis[['d2haz']] <- d2haz
  lis[['dcumhaz']] <- dcumhaz
  lis[['d2cumhaz']] <- d2cumhaz

  return(lis)

}

input_se_functions <- function( groups, time, status, X, beta_hat, haz, cumhaz, dhaz, dcumhaz, d2haz, d2cumhaz )
{
  H <- length( unique( groups ) )
  ncovs <- length( beta_hat )


  lamroutlam0 <- matrix( nrow = H, ncol = ncovs )
  sumhaz <- rep( 0, H )
  dsumhaz <- array( dim = c( H, ncovs ) )
  d2sumhaz <- lamsum <- array( dim = c( H, ncovs, ncovs ) )
  deltaijcovr <- matrix( 0, nrow = H, ncol = ncovs )
  dsumhazred <- matrix( 0, nrow = H, ncol = ncovs )


  for( j in 1:H )
  {
    hospj <- groups==j
    ebz <- exp( as.matrix( X[hospj,,drop=FALSE] ) %*% beta_hat )
    ti <- match(time[hospj], cumhaz$time)
    lam0t <- cumhaz$hazard[ti]
    sumhaz[j] <- sum(lam0t * ebz)


    for (r in 1:ncovs){
      covr <- X[hospj,r]
      lamrt <- dcumhaz[ti, r]
      #dsumhaz is the second part of III term
      dsumhaz[j,r] <- sum(ebz*(covr*lam0t + lamrt))
      #new lines for III term
      lamroutlam0tmp <- status[hospj]*( dhaz[ti,r]/haz[ti] + covr )
      lamroutlam0[j,r] <- sum(lamroutlam0tmp[haz[ti]>0])

      deltaijcovrtmp <- status[hospj]*covr
      deltaijcovr[j,r] <- sum(deltaijcovrtmp[haz[ti]>0])
      dsumhazred[j,r] <- sum(ebz*(covr*lam0t))

      for (s in 1:ncovs){
        covs <- X[hospj,s]
        lamst <- dcumhaz[ti, s]
        lamrst <- d2cumhaz[ti, r, s]
        d2sumhaz[j,r,s] <- sum(ebz*(covr*covs*lam0t + covr*lamst + covs*lamrt + lamrst))
        lamtmp <- status[hospj] * (d2haz[ti,r,s]*haz[ti] - dhaz[ti,r]*dhaz[ti,s])/haz[ti]^2
        lamsum[j,r,s] <- sum(lamtmp[haz[ti]>0]) # exclude zero hazard

      }

    }
  }

  lis <- list()
  lis[['lamroutlam0']] <- lamroutlam0
  lis[['lamsum']] <- lamsum
  lis[['deltaijcovr']] <- deltaijcovr
  lis[['sumhaz']] <- sumhaz
  lis[['dsumhazred']] <- dsumhazred
  lis[['dsumhaz']] <- dsumhaz
  lis[['d2sumhaz']] <- d2sumhaz

  return(lis)
}

log_likelihood <- function( groups, time, status, X, p, w, alpha, beta_hat, haz, cumhaz)
{
  H <- dim(alpha)[1]
  K <- dim(alpha)[2]
  ncovs <- ncol(X)

  llik = 0

  for( j in 1:H )
  {
    hospj <- groups==j
    hospj_lam <- groups==j & status==1
    ti <- match(time[hospj], cumhaz$time)
    ti_lam <- match(time[hospj_lam], cumhaz$time)

    for( k in 1:K )
    {
      llik_LAM <- ifelse( ncovs > 0,
                          sum( cumhaz$hazard[ti]*w[k]*exp( as.matrix( X[hospj,,drop=FALSE] ) %*% beta_hat ) ),
                          sum( cumhaz$hazard[ti]*w[k] ) )

      #w_off[hospj][1] because w_off[hospj] are all equal
      llik_lam <- ifelse( ncovs > 0,
                          sum( log( haz[ti_lam]*w[k]*exp( as.matrix( X[hospj_lam,,drop=FALSE] ) %*% beta_hat ) ) ),
                          sum( log( haz[ti_lam]*w[k] ) ) )

      llik_p <- ifelse( K==1, 1, p[ k ] )

      llik <- llik + alpha[j,k]*(llik_lam-llik_LAM) + alpha[j,k]*log( llik_p )

    }
  }
  return(llik)
}

