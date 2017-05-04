#' Estimating a Cox model with a shared nonparametric frailty
#'
#' This function estimates the parameter of a semiparametric Cox model with a nonparametric frailty term.
#'
#' @param data a data.frame in which there are status, time, family and covariates variables
#' @param time name of the variable which indicates the status
#' @param status name of the variable which indicates the status
#' @param family name of the variable by which the data are divided
#' @param covariates vector of covariates' names
#' @param population initial number of latent populations
#' @param eps_conv tolerance
#'
#' @return This fuction returns a list of elements, whose length is equal to the initial number of populations + 4. The last four elements of the list show the best model according to loglikelihood, BIC and AIC, while the fourth term represents the minimum number of latent population for which each population has at least one member.
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
#' data <- simulWeibDiscreteFrailCovNPbaseInv( N, S, Lambda_0_inv, rho, beta, p, w_values)
#'
#' test_res <- npdf_cox( data, time = 'time', status = 'status', family = 'family', covariates = c('x'), population = 4, eps_conv=10^-4)
#' best_model_llik <- test_res[[5]] #population+1
#' best_model_BIC <- test_res[[6]] #population+2
#' best_model_AIC <- test_res[[7]] #population+3
#' best_model_K <- test_res[[8]] #population+4
#'
#' test_res[[best_model_BIC]]



npdf_cox <- function( data, time, status, family, covariates, population = 2, eps_conv=10^-4)
{
  H <- length( unique( data[[family]] ) ) #total number of groups
  K <- population #number of latent populations
  dat_ord <- data[ order( data[[family]], data[[time]], -data[[status]] ), ]
  nt <- length(unique(data[[time]]))
  ncovs <- length(covariates)
  cov_read <- paste(covariates, sep=',')
  cumhaz <- as.data.frame( cbind( hazard=rep( 0, nt ),
                                   time = sort( unique( data[[time]] ) )) )

  #### computing SE
  YY <- haz <- rep( 0, nt )

  ##sort the cluster id (there is a smarter way)
  for( i in 1:length(unique( dat_ord[[family]] ) ) )
  {
    value <- unique( dat_ord[[family]] )[i]
    dat_ord[[family]][ dat_ord[[family]] == value ] <- i
  }

  D <- table( factor( dat_ord[[family]] )[ dat_ord[[status]] == 1 ] )  # number of events in each group j
  n_H <- as.vector( table( dat_ord[[family]] ) ) #number of patients in each hospital

  risk_index <- matrix( 0, nrow = dim(dat_ord)[1], ncol = nt ) # matrix of at risk patients: nrows x nt
  N = rep( 0, nt )

  time_list <- sapply( 1:nt, function(x) !is.na(match(dat_ord[[time]], cumhaz$time[x]))) #matrix dim(dat_ord)[1] x nt, time_list[i,j] = T if the i-th row in the data is the j-th time recorded
  N <- sapply( 1:dim(time_list)[2], function(x) sum(dat_ord[[status]][time_list[,x]])) #number of events happenend at each time. N[j] is the number of eents happenend at time j


  for( l in 1:nt )
  {
    risk_index[ which( dat_ord[[time]] >= cumhaz$time[ l ]), l ] <-  1 #nrow(dat_ord) x nt
  }

  lis = list()

  while( K >= 1 )
  {
    count <- 0
    not_proper_K <- 0
    eps <- 10^5

    p <- rep( 1/K, K )
    w <- seq( 1, 5, length.out = K )

    if( K==1 ) w <- 1

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

      # estimating the betas
      beta_hat <- temp_model$coef

      # computing baseline hazard and cumulative baseline hazard
      YY = as.numeric( (w_off*exp( beta_hat %*% t(dat_ord[,cov_read]) )) %*% risk_index ) # nt x 1 contribute from at risk patient
      haz = N/YY
      cumhaz$hazard = cumsum( haz )

      # saving p
      p_old <- p


      #Expectation step

      for( j in 1:H )
      {
        hospj <- dat_ord[[family]]==j
        ebz <- exp( as.matrix( dat_ord[hospj,cov_read] ) %*% beta_hat )
        ti <- match(dat_ord[[time]][hospj], cumhaz$time)
        lam0t <- cumhaz$hazard[ti]


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

      # vector of frailty, dim = Kx1
      w <- (D %*% alpha)/(E_formula_part %*% alpha)

      # standardization of w
      w <- w/min(w)
      #w[which.min(w)] <- 1

      # vector of frailty values, dim = nx1
      w_off <- ( w[ belonging ][dat_ord[[family]]] )

      # computation of the distance from the previous estimate of p
      eps <- max( abs(p-p_old) )

      # augmenting count
      count <- count + 1
    }

    # computing 1st and 2nd hazard derivatives
    dhazards <- harzard_derivatives( dat_ord, risk_index, w_off, beta_hat, cov_read, N, YY )
    dhaz <- dhazards$dhaz
    d2haz <- dhazards$d2haz
    dcumhaz <- dhazards$dcumhaz
    d2cumhaz <- dhazards$d2cumhaz

    #computing the input for the se functions
    input_se <- input_se_functions( family, time, status, dat_ord, cov_read, beta_hat, haz, cumhaz, dhaz, dcumhaz, d2haz, d2cumhaz )
    lamroutlam0 <- input_se$lamroutlam0
    lamsum <- input_se$lamsum
    deltaijcovr <- input_se$deltaijcovr
    sumhaz <- input_se$sumhaz
    dsumhazred <- input_se$dsumhazred
    dsumhaz <- input_se$dsumhaz
    d2sumhaz <- input_se$d2sumhaz


    #function to be numerically optimized
    sefundeponwbeta = function( x )
    {
      YYfun <- as.numeric( (w_off*exp( x[(2*K):(2*K-1+ncovs)] %*% t(dat_ord[,cov_read]) )) %*% risk_index ) # nt x 1 contribute from at risk patient
      hazfun <- N/YYfun
      cumhazfun <- cumsum(hazfun)

      lfull <- rep(0,H)
      for(j in 1:H)
      {
        hospj <- dat_ord[[family]]==j
        ebz <- exp( as.matrix( dat_ord[hospj,cov_read] ) %*% x[(2*K):(2*K-1+ncovs)] )
        ti <- match(dat_ord[[time]][hospj], cumhaz$time)
        lam0t <- cumhazfun[ti]
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

      #Louis method for standard errors
      EStS <- Ixy( ncovs, alpha, w, p, D, sumhaz, lamroutlam0, dsumhaz)
      EB <- Ix( ncovs, alpha, w, p, D, d2sumhaz, lamsum, dsumhaz)
      ESES <- Iy( ncovs, alpha, w, p, D, sumhaz, lamroutlam0, dsumhaz)
      stderr_covariance <- information_matrix( expectedS = EStS, expectedB = EB, expectedSexpectedS = ESES )

      #Exact second derivativefor standard errors
      SecDer <- exact_sec_der( ncovs, w, p, D, sumhaz, dsumhaz, d2sumhaz,lamsum )
      InfoFischer <- solve( -SecDer )
    }


    # computing loglikelihood, AIC and BIC
    llik <- log_likelihood( family, time, status, dat_ord, cov_read, p, w, alpha, beta_hat, haz, cumhaz)

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
  #Using the right K as the top value for evaluating the lowest BIC and AIC
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
