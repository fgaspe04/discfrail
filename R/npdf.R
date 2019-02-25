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

#' Cox model for grouped survival data with nonparametric discrete shared frailties
#'
#' This function fits a Cox proportional hazards model to grouped survival data, where the shared group-specific frailties have a nonparametric discrete distribution.  An EM algorithm is used to maximise the marginal partial likelihood.
#'
#' @param formula A formula expression in conventional R linear modelling
#' syntax. The response must be a survival time constructed by the
#' \code{\link{Surv}} function from the \pkg{survival} package, and
#' any covariates are given on the right-hand
#' side.  For example,
#'
#' \code{Surv(time, dead) ~ age + sex}
#'
#' Only \code{Surv} objects of \code{type="right"} are supported, corresponding to right-censored observations.
#'
#' @param groups name of the variable which indicates the group in which each individual belongs (e.g. the hospital that the individual is treated in).  This can be integer, factor or character.  The name should be unquoted.
#'
#' @param data A data frame in which to find variables supplied in
#' \code{formula}.  If not given, the variables should be in the working
#' environment.
#'
#' @param K initial number of latent populations, or clusters of groups which have the same discrete frailty.
#'
#' @param estK If \code{TRUE} (the default) then
#' multiple models are fitted with number of latent groups ranging from 1 to \code{K}.  The "best fitting" model according to the criterion specified in \code{criterion} is then highlighted when printing the object returned by this function.
#'
#' If \code{FALSE} then the number of latent populations is fixed at \code{K}.
#'
#' @param criterion Criterion used to choose the best-fitting model to highlight when \code{estK} is \code{TRUE}.
#'
#' \code{"Laird"} for the Laird criterion. Running from K latent populations to 1 latent population, this criterion selects the maximum number of latent populations that are non empty as the best K.
#'
#' \code{"AIC"} for Akaike's information criterion.
#'
#' \code{"BIC"} for the Bayesian information criterion (the default).
#'
#' @param eps_conv convergence tolerance for the EM algorithm.
#'
#' @param se_method Method or methods used to compute the standard errors.  A character vector containing one or more of the following:
#'
#' \code{"louis"} The method of Louis (1982) based on an approximation to the information matrix.
#'
#' \code{"exact"} In this method the standard errors are computed directly from the observed information matrix obtained by analytic differentiation.
#'
#' \code{"numeric"} This method uses numerical differentiation to approximate the information matrix, and is substantially slower.
#'
#' By default this is \code{c("louis","exact")} because these two methods are equally fast. So that SEs from both these two methods are calculated and presented.   Set \code{se_method=NULL} to compute no standard errors.
#'
#' @return If \code{estK=FALSE} this returns a list of class \code{npdf} which includes information about the model fit, including estimates and standard errors.
#'
#' If \code{estK=TRUE} this returns a list of class \code{npdflist}.  This has an element \code{models} that contains a list of length \code{K}, with one component of class \code{npdf} for each fitted model.
#'
#'\code{comparison} is a matrix composed of \code{K} rows and 5 columns (\code{K}, \code{K_fitted}, \code{llik}, \code{AIC}, \code{BIC}). \code{K_fitted} is the number of estimated latent populations, which can be equal to or less than \code{K}. \code{llik} stands for log-likelihood, \code{AIC} for Akaike Information Criterion and \code{BIC} for Bayesian Information Criterion.
#'
#'\code{Kopt} is optimal model under each criterion.
#'
#'\code{criterion} is the preferred criterion.
#'
#' In either case, the data frame used for the fit (the "model frame") is appended as a component \code{mf}.
#'
#'
#' @references
#' Gasperoni, F., Ieva, F., Paganoni, A.M., Jackson, C. and Sharples, L. (2018). Nonparametric frailty Cox models for hierarchical time-to-event data. \emph{Biostatistics}.
#'
#' Laird, N. (1978). Nonparametric maximum likelihood estimation of a mixing distribution. \emph{Journal of the American Statistical Association}, 73(364), 805–811.
#'
#' Louis, T. A. (1982). Finding the observed information matrix when using the EM algorithm. \emph{Journal of the Royal Statistical Society. Series B}, 44(2), 226–233.
#'
#' @import survival
#' @importFrom stats model.extract model.matrix rnorm rpois runif
#' @importFrom utils tail
#' @importFrom numDeriv genD
#' @importFrom Matrix forceSymmetric
#'
#' @export
#'
#' @examples
#' test <- npdf_cox( Surv(time, status) ~ x, groups=family, data=weibdata2030, K = 4, eps_conv=10^-4)
#' test    # optimal model (by all criteria) has 2 latent populations
#' test$models[[1]] # examine alternative model with 1 latent population
#'
npdf_cox <- function(formula, groups, data,
                     K=2, estK=TRUE, criterion="BIC",
                     eps_conv=10^-4, se_method=c("louis","exact")){
  call <- match.call()
  indx <- match(c("formula", "groups", "data"), names(call), nomatch = 0)
  if (indx[1] == 0)
    stop("A \"formula\" argument is required")
  if (indx[2] == 0)
    stop("A \"groups\" argument is required")
  temp <- call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  temp[["formula"]] <- formula
  if (missing(data)) temp[["data"]] <- environment(formula)
  mf <- eval(temp, parent.frame())

  Y <- model.extract(mf, "response")
  time <- Y[,"time"]
  status <- Y[,"status"]
  groups <- model.extract(mf, "groups")
  mf <- mf[order(groups, time, -status),,drop=FALSE] # is the ordering necessary?
  Y <- model.extract(mf, "response")
  time <- Y[,"time"]
  status <- Y[,"status"]
  groups <- model.extract(mf, "groups")

  mm <- model.matrix(formula, mf)
  X <- mm[,-1,drop=FALSE]

  H <- length( unique( groups ) ) #total number of groups
  nt <- length(unique(time))

  nobs <- nrow(mf)
  ncovs <- ncol(X)
  cumhaz <- as.data.frame( cbind( hazard=rep( 0, nt ),
                                 time = sort( unique( time ) )) )

  #### computing SE
  YY <- haz <- rep( 0, nt )

  ## convert the cluster id to numeric
  groups <- match(groups, unique(groups))

  D <- table( factor( groups )[ status == 1 ] )  # number of events in each group j
  risk_index <- matrix( 0, nrow = nobs, ncol = nt ) # matrix of at risk patients: nrows x nt

  time_list <- sapply( 1:nt, function(x) !is.na(match(time, cumhaz$time[x]))) # matrix dim(dat_ord)[1] x nt, time_list[i,j] = T if the i-th row in the data is the j-th time recorded
  N <- sapply( 1:dim(time_list)[2], function(x) sum(status[time_list[,x]])) # N[j] is the number of events which happenened at time j

  for( l in 1:nt )
  {
    risk_index[ which( time >= cumhaz$time[ l ]), l ] <-  1 #nrow(dat_ord) x nt
  }

  if (estK) {
      models <- vector(K, mode="list")
      k <- K
      while( k >= 1 ) {
          models[[k]] <- npdf_core(formula=formula, data=data, K=k,
                                   time=time, status=status, groups=groups,
                                   X=X, ncovs=ncovs, N=N, H=H, cumhaz=cumhaz,
                                   D=D, risk_index = risk_index, eps_conv=eps_conv,
                                   se_method = se_method)
          k <- k-1
      }
      ## TODO disagree should truncate AIC selection above at Laird best fit K

      ## Matrix of fit statistics for each model
      comparison <- t(sapply(models,
                           function(x) unlist(x[c("K_fitted","llik","AIC", "BIC")])
                           ))
      comparison <- as.data.frame(cbind(K=1:K, comparison))
      Kopt <- c("Laird" = max(which(comparison$K == comparison$K_fitted)),
                "AIC" = which.min(comparison$AIC),
                "BIC" = which.min(comparison$BIC))
      res <- list(models=models, comparison=comparison, Kopt=Kopt, criterion=criterion)
      class(res) <- "npdflist"
  }  else  {
      res <- npdf_core(formula=formula, data=data, K=K,
                       time=time, status=status, groups=groups,
                       X=X, ncovs=ncovs, N=N, H=H, cumhaz=cumhaz,
                       D=D, risk_index = risk_index, eps_conv=eps_conv,
                       se_method = se_method)
  }
  res$call <- call
  res$mf <- mf
  res
}


npdf_core <- function(formula,
                      data=NULL,
                      K=2,
                      time,
                      status,
                      groups,
                      X,
                      ncovs,
                      N,
                      H,
                      cumhaz,
                      D,
                      risk_index,
                      eps_conv=10^-4,
                      se_method
                      ){

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

    w_off <- w_to_rep[ groups ]

    numerator <- rep( 0, K )
    E_formula_part <- rep( 0, H )
    alpha <- E_formula <- matrix( 0, nrow = H, ncol = K)

    while(eps > eps_conv & count < 200 )
    {

      form_off <- formula( paste(deparse(formula), "offset( log( w_off ) )", sep ='+'), sep = "~")
      temp_model <- coxph( formula=form_off, data=data, method = "breslow" )

      # estimating the betas
      beta_hat <- temp_model$coef

      # computing baseline hazard and cumulative baseline hazard
      YY = as.numeric( (w_off*exp( beta_hat %*% t(X) )) %*% risk_index ) # nt x 1 contribute from at risk patient
      haz = N/YY
      cumhaz$hazard = cumsum( haz )

      # saving p
      p_old <- p


      #Expectation step

      for( j in 1:H )
      {
        hospj <- groups==j
        ebz <- exp( as.matrix( X[hospj,] ) %*% beta_hat )
        ti <- match(time[hospj], cumhaz$time)
        lam0t <- cumhaz$hazard[ti]

        E_formula_part[j] <- ifelse( ncovs > 0,
                                     sum( lam0t*ebz ),
                                     sum( lam0t ) )

        for( k in 1:K )
        {
          E_formula[j,k] <- ifelse( ncovs > 0,
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
      w_off <- ( w[ belonging ][groups] )

      # computation of the distance from the previous estimate of p
      eps <- max( abs(p-p_old) )

      # augmenting count
      count <- count + 1
    }


    # computing 1st and 2nd hazard derivatives
    dhazards <- hazard_derivatives( X, risk_index, w_off, beta_hat, N, YY )
    dhaz <- dhazards$dhaz
    d2haz <- dhazards$d2haz
    dcumhaz <- dhazards$dcumhaz
    d2cumhaz <- dhazards$d2cumhaz

    #computing the input for the se functions
    input_se <- input_se_functions( groups, time, status, X, beta_hat, haz, cumhaz, dhaz, dcumhaz, d2haz, d2cumhaz )
    lamroutlam0 <- input_se$lamroutlam0
    lamsum <- input_se$lamsum
    deltaijcovr <- input_se$deltaijcovr
    sumhaz <- input_se$sumhaz
    dsumhazred <- input_se$dsumhazred
    dsumhaz <- input_se$dsumhaz
    d2sumhaz <- input_se$d2sumhaz


    #function to be numerically optimized
    sefundeponwbeta = function( x ) #length(x) = K-1 + K + ncovs
    {
      YYfun <- as.numeric( (w_off*exp( x[(2*K):(2*K-1+ncovs)] %*% t(X) )) %*% risk_index ) # nt x 1 contribute from at risk patient #x[(2*K):(2*K-1+ncovs)] are betas
      hazfun <- N/YYfun
      cumhazfun <- cumsum(hazfun)

      lfull <- rep(0,H)
      for(j in 1:H)
      {
        hospj <- groups==j
        ebz <- exp( as.matrix( X[hospj,] ) %*% x[(2*K):(2*K-1+ncovs)] )
        ti <- match(time[hospj], cumhaz$time)
        lam0t <- cumhazfun[ti]
        haz0t <- hazfun[ti]

        firsttermtmp <- status[hospj]*log(haz0t*ebz)
        firstterm <- sum(firsttermtmp[haz0t>0,])

        secondtermtmp <- -sum(lam0t * ebz) * x[(K):(2*K-1)] #x[(K):(2*K-1)] are w_1,..,w_K
        wpowerD <- x[(K):(2*K-1)]^D[j]
        secondterm <- log(sum(exp(secondtermtmp)*wpowerD*c(x[1:(K-1)],1-sum(x[1:(K-1)])))) #x[1:(K-1)] are pi_1,..,\pi_{K-1}

        lfull[j] <- firstterm+secondterm
      }
      return( sum( lfull ) )
    }

    # computing loglikelihood, AIC and BIC
    llik <- log_likelihood( groups, time, status, X, p, w, alpha, beta_hat, haz, cumhaz)

    tot_param <- K-1 + K + ncovs
    BIC <- -2*llik + tot_param * log( length( which(status == 1 ) ) )
    AIC <- -2*llik + 2*tot_param

    K_fitted <- length(unique(belonging))
    if( K_fitted < K )
    {
      not_proper_K <- 1
      stderr_covariance <- NULL
      infoRich <- NULL
      InfoFischer <- NULL
    }
    names(p) <- paste0("p",1:K)
    names(w) <- paste0("w",1:K)

    res <- list(
        K = K,
        too_much_K = not_proper_K,
        K_fitted = K_fitted,
        p = p,
        w = w,
        beta = beta_hat,
        belonging = belonging,
        llik = llik,
        BIC = BIC,
        AIC = AIC,
        alpha = alpha
    )

    if( K_fitted == K )
    {
        if ("numeric" %in% se_method){
            ##Numerical computation of standard errors
            if( K != 1 )
            {
                deriv <- numDeriv::genD( sefundeponwbeta, c( p[1:(K-1)], w, beta_hat ), method = "Richardson" )
                info <- matrix(0, 2*K-1+ncovs, 2*K-1+ncovs)
                info[upper.tri(info, diag=T)] <- deriv$D[(2*K+ncovs):(2*K-1+ncovs+(2*K+ncovs-1)*(2*K+ncovs)/2)]
                info <- Matrix::forceSymmetric(info)
                infoNumeric <- solve(-info)
                res$varcovNumeric = infoNumeric
                res$seNumeric = NULL
                res$seNumeric[ 1:( K - 1 ) ] = sqrt(diag(infoNumeric))[ 1:( K - 1 ) ]
                res$seNumeric[ K ] = sqrt( sum( infoNumeric[ ( 1:( K - 1 ) ), ( 1:( K - 1 ) ) ] ) )
                res$seNumeric[ ( K + 1 ):( 2*K - 1 ) ] = sqrt( (w[-1]/w[1])^2*((sqrt(diag(infoNumeric))[ (K + 1):(2*K - 1) ]/w[-1])^2+( sqrt(diag(infoNumeric))[ K ]/w[1])^2-2*infoNumeric[ (K + 1):(2*K - 1), K ]/(w[-1]*w[1])) )
                res$seNumeric[ (2*K):( 2*K + ncovs - 1 ) ] = sqrt(diag(infoNumeric))[ (2*K):( 2*K + ncovs - 1 ) ]
            }else{
                deriv <- numDeriv::genD( sefundeponwbeta, c( w, beta_hat), method = "Richardson" )
                info <- matrix(0, 2*K-1+ncovs, 2*K-1+ncovs)
                info[upper.tri(info, diag=T)] <- deriv$D[(2*K+ncovs):(2*K-1+ncovs+(2*K+ncovs-1)*(2*K+ncovs)/2)]
                info <- Matrix::forceSymmetric(info)
                infoNumeric <- solve(-info)
                res$varcovNumeric = infoNumeric
                res$seNumeric = NULL
                res$seNumeric[ 1: (ncovs) ] = sqrt(diag(infoNumeric))[ 2:( ncovs + 1 ) ]
            }

        }

        if ("louis" %in% se_method){
            ##Louis method for standard errors
            EStS <- Ixy( ncovs, alpha, w, p, D, sumhaz, lamroutlam0, dsumhaz)
            EB <- Ix( ncovs, alpha, w, p, D, d2sumhaz, lamsum, dsumhaz)
            ESES <- Iy( ncovs, alpha, w, p, D, sumhaz, lamroutlam0, dsumhaz)
            stderr_covariance <- information_matrix( expectedS = EStS, expectedB = EB, expectedSexpectedS = ESES )
            res$varcovLouis = stderr_covariance[[1]]
            res$seLouis = NULL
            if( K != 1 ){
              res$seLouis = NULL
              res$seLouis[ 1:( K - 1 ) ] = stderr_covariance[[2]][ 1:( K - 1 ) ]
              res$seLouis[ K ] = sqrt( sum( stderr_covariance[[1]][ ( 1:( K - 1 ) ), ( 1:( K - 1 ) ) ] ) )
              res$seLouis[ ( K + 1 ):( 2*K - 1 ) ] = sqrt( (w[-1]/w[1])^2*((stderr_covariance[[2]][ (K + 1):(2*K - 1) ]/w[-1])^2+( stderr_covariance[[2]][ K ]/w[1])^2-2*stderr_covariance[[1]][ (K + 1):(2*K - 1), K ]/(w[-1]*w[1])) )
              res$seLouis[ (2*K):( 2*K + ncovs - 1 ) ] = stderr_covariance[[2]][ (2*K):( 2*K + ncovs - 1 ) ]
            }else{
              res$seLouis[ 1:( ncovs ) ] = stderr_covariance[[2]][ 2:( ncovs + 1 ) ]
            }
        }
        if ("exact" %in% se_method){
            ##Exact second derivativefor standard errors
            SecDer <- exact_sec_der( ncovs, w, p, D, sumhaz, dsumhaz, d2sumhaz,lamsum )
            InfoFisher <- solve( -SecDer )
            res$varcovExact = InfoFisher
            res$seExact = NULL
            if( K != 1 ){
              res$seExact = NULL
              res$seExact[ 1:( K - 1 ) ] = sqrt(diag(InfoFisher))[ 1:( K - 1 ) ]
              res$seExact[ ( K ) ] = sqrt( sum( InfoFisher[ ( 1:( K - 1 ) ), ( 1:( K - 1 ) ) ] ) )
              res$seExact[ ( K + 1 ):( 2*K - 1 ) ] = sqrt( (w[-1]/w[1])^2*((sqrt(diag(InfoFisher))[ (K + 1):(2*K - 1) ]/w[-1])^2+( sqrt(diag(InfoFisher))[ K ]/w[1])^2-2*InfoFisher[ (K + 1):(2*K - 1), K ]/(w[-1]*w[1])) )
              res$seExact[ (2*K):( 2*K + ncovs - 1 ) ] = sqrt(diag(InfoFisher))[ (2*K):( 2*K + ncovs - 1 ) ]
            }else{
              res$seExact[ 1:( ncovs ) ] = sqrt(diag(InfoFisher))[ 2:( ncovs + 1 ) ]
            }
        }
    }


    class(res) <- "npdf"
    res
}
