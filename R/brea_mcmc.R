
brea_mcmc <- function(x, y, priors = NULL, S = 1000, B = 100, n = NULL,
                      K = NULL, store_re = FALSE, block_mh = TRUE) {
  
  
  # verify and clean user input ------------------------------------------------

  
  # verify and set x and K:
  xK <- check_set_predictors(x,K)
  x <- xK[[1]]
  K <- xK[[2]]
  
  # set the number of predictors:
  M <- ncol(x)
  
  # set the number of observations:
  N <- nrow(x)
  
  # verify and set y and n:
  yn <- check_set_outcome(y,n,N)
  y <- yn[[1]]
  n <- yn[[2]]
  
  # set the number of event types (i.e. competing risks):
  R <- ncol(y)

  # verify and set prior specifications:
  priors <- check_set_priors(priors,M,R)
  
  # verify number of iterations to store S is a positive whole number:
  if (!(is.numeric(S) && length(S)==1 && all_whole(S) && S > .5)) {
    stop("the number of iterations to store S must be a positive whole number")
  }
  
  # convert S to integer type if necessary:
  if(!is.integer(S)) S <- as.integer(round(S))
  
  # verify number of burn-in iterations B is a nonnegative whole number:
  if (!(is.numeric(B) && length(B)==1 && all_whole(B) && B > -.5)) {
    stop("the number of burn-in iterations B must be a nonnegative integer")
  }
  
  # convert B to integer type if necessary:
  if(!is.integer(B)) B <- as.integer(round(B))

  
  # initialize parameters, parameter storage, and acceptance counters ----------

  
  # initialize current values of intercepts to MLE's with covariates omitted:
  b_0 <- log(colSums(y)/(sum(n)-sum(y)))
  
  # intercept storage:
  b_0_s <- matrix(0,R,S)
  
  # current values, storage and acceptance counters for linear predictor params:
  b_m <- vector("list",M)
  b_m_s <- vector("list",M)
  b_m_a <- vector("list",M)
  
  for (m in seq_len(M)) {
    # current parameter values:
    b_m[[m]] <- matrix(0,R,K[m])
    # stored parameter values (saving random effects only if requested):
    if (priors[[m]]$type %in% c("cat","gmrf") ||
        (priors[[m]]$type == "re" && store_re)) {
      b_m_s[[m]] <- array(0,c(R,K[m],S))
    }
    # acceptance counters:
    if (priors[[m]]$type == "gmrf" && block_mh) b_m_a[[m]] <- 0L
    else b_m_a[[m]] <- rep(0L,K[m])
  }
  
  # current values and storage for variance/precision parameters:
  
  # std. dev. of the random walk increment or covariance of the random effects:
  s_m <- vector("list",M)
  s_m_s <- vector("list",M)
  # precision of the random walk increment or precision of the random effects:
  prec_m <- vector("list",M)

  for (m in seq_len(M)) {
    
    # GMRF case:
    if (priors[[m]]$type == "gmrf") {
      # initialize increment std. dev. to 1:
      s_m[[m]] <- rep(1,R)
      s_m_s[[m]] <- matrix(0,R,S)
      prec_m[[m]] <- 1/s_m[[m]]^2
      
    # random effects case:
    } else if (priors[[m]]$type == "re") {
      # initialize random effects covariance to identity matrix:
      s_m[[m]] <- diag(R)
      s_m_s[[m]] <- array(0,c(R,R,S))
      prec_m[[m]] <- solve(s_m[[m]])
    }
    
  }
  
  # initialize linear predictors and negative log partition function:
  
  # current values of the linear predictors:
  eta <- matrix(0,R,N)
  eta <- eta + b_0
  for (m in seq_len(M)) eta <- eta + b_m[[m]][,x[,m]]
  
  # current values of the negative log partition functions:
  nlpartfun <- -log1p(colSums(exp(eta)))
  
  
  # create index subsets and y-sums for single-category RWM steps --------------
  
  
  # list whose m^th element is the list whose k^th element is the vector of
  # indices of observations with value k for predictor m:
  subindex <- vector("list",M)
  
  # list whose m^th element is the R x K[m] matrix whose (r,k) entry is the
  # count of events of type r occurring at observations with value k for
  # predictor m:
  y_sum <- vector("list",M)
  
  # loop over predictors m:
  for (m in seq_len(M)) {
    subindex[[m]] <- vector("list",K[m])
    y_sum[[m]] <- matrix(0,R,K[m])
    for (k in seq_len(K[m])) {
      subindex[[m]][[k]] <- which(x[,m]==k)
      y_sum[[m]][,k] <- colSums(y[subindex[[m]][[k]],,drop=FALSE])
    }
  }
  
  
  # run MCMC -------------------------------------------------------------------

  
  # loop over scans s from both burn-in and stored iterations:  
  for (s in seq_len(S+B)) {
  
    # loop over predictors m:  
    for (m in seq_len(M)) {
      
      
      # categorical case RWM updates -------------------------------------------
      
      
      if (priors[[m]]$type == "cat") {
        
        for (k in seq_len(K[m])) {
          
          # proposed step in the random walk:
          prop_step <- rnorm(R,0,
                             2.4/sqrt(R*(priors[[m]]$tau + y_sum[[m]][,k])))
          
          # proposed new parameter values:
          b_star <- b_m[[m]][,k] + prop_step
          
          # conditional prior means:
          cm <- rowMeans(b_m[[m]][,-k,drop=FALSE])
          
          # calculate log prior ratio:
          lpr <- sum(dnorm(b_star,cm,priors[[m]]$sd,TRUE)
                     - dnorm(b_m[[m]][,k],cm,priors[[m]]$sd,TRUE))
          
          # index subset of observations with value k for variable m:
          subi <- subindex[[m]][[k]]
          
          # if there are no observations with value k for variable m, then
          # accept the proposal based on the log prior ratio only:
          if (length(subi) == 0L) {
            
            if (log(runif(1)) < lpr) {
              b_m[[m]][,k] <- b_star
              if (s > B) b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
            }

          # otherwise, calculate the log likelihood ratio then accept/reject:
          } else {
            
            # linear predictors under proposed values:
            eta_star <- eta[,subi] + prop_step
            
            # negative log partition function under proposed values:
            nlpartfun_star <- -log1p(.colSums(exp(eta_star),R,length(subi)))
            
            # log likelihood ratio:
            llr <- (crossprod(y_sum[[m]][,k],prop_step) + 
                      crossprod(n[subi],nlpartfun_star-nlpartfun[subi]))
            
            # accept the proposal with the appropriate probability:
            if (log(runif(1)) < (llr + lpr)) {
              b_m[[m]][,k] <- b_star
              eta[,subi] <- eta_star
              nlpartfun[subi] <- nlpartfun_star
              if (s > B) b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
            }
            
          }
          
        } # end k loop
        
        # sweep parameter vector averages onto intercepts:
        rms <- rowMeans(b_m[[m]])
        b_m[[m]] <- b_m[[m]] - rms
        b_0 <- b_0 + rms
        
      
      # GMRF case updates ------------------------------------------------------
      
        
      } else if (priors[[m]]$type == "gmrf") {
        
        # RWM case -------------------------------------------------------------

        # perform RWM updates pre-burn-in and if block M-H updates not used:
        if (s <= B || !block_mh) {
          
          for (k in seq_len(K[m])) {
          
            # approximate full conditional precision:
            fc_prec <- ifelse(k==1L || k==K[m],1,2)*prec_m[[m]] + y_sum[[m]][,k]
            
            # proposed step in the random walk:
            prop_step <- rnorm(R,0,2.4/sqrt(R*fc_prec))
            
            # proposed new parameter values:
            b_star <- b_m[[m]][,k] + prop_step
            
            # log prior ratio:
            if (k==1L) {
              lpr <- sum(dnorm(b_star,b_m[[m]][,2L],s_m[[m]],TRUE) -
                           dnorm(b_m[[m]][,k],b_m[[m]][,2L],s_m[[m]],TRUE))
            } else if (k==K[m]) {
              lpr <- sum(dnorm(b_star,b_m[[m]][,k-1L],s_m[[m]],TRUE) -
                           dnorm(b_m[[m]][,k],b_m[[m]][,k-1L],s_m[[m]],TRUE))
            } else {
              b_near <- (b_m[[m]][,k-1L] + b_m[[m]][,k+1L])/2
              lpr <- sum(dnorm(b_star,b_near,s_m[[m]]/sqrt(2),TRUE) -
                           dnorm(b_m[[m]][,k],b_near,s_m[[m]]/sqrt(2),TRUE))
            }
            
            # index subset of observations with value k for variable m:
            subi <- subindex[[m]][[k]]
            
            # if there are no observations with value k for variable m, then
            # accept the proposal based on the log prior ratio only:
            if (length(subi) == 0L) {
              
              if (log(runif(1)) < lpr) {
                b_m[[m]][,k] <- b_star
                if (s > B) b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
              }
              
            # otherwise, calculate the log likelihood ratio then accept/reject:
            } else {
              
              # linear predictors under proposed values:
              eta_star <- eta[,subi] + prop_step
              
              # negative log partition function under proposed values:
              nlpartfun_star <- -log1p(.colSums(exp(eta_star),R,length(subi)))
              
              # log likelihood ratio:
              llr <- (crossprod(y_sum[[m]][,k],prop_step) + 
                        crossprod(n[subi],nlpartfun_star-nlpartfun[subi]))
              
              # accept the proposal with the appropriate probability:
              if (log(runif(1)) < (llr + lpr)) {
                b_m[[m]][,k] <- b_star
                eta[,subi] <- eta_star
                nlpartfun[subi] <- nlpartfun_star
                if (s > B) b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
              }
              
            }
          
          } # end k loop

        # block Metropolis-Hastings case ----------------------------------------

        } else {
          
          # step 1: sample from approximate full conditional -------------------

          # exponentials of linear predictors (R x N):
          e_eta <- exp(eta)
          
          # discrete hazards for each observation and risk (N x R):
          lm <- t(e_eta)/(1+.colSums(e_eta,R,N))
          
          # discrete hazards adjusted by replication count (N x R):
          lm_n <- n*lm
          
          # hazard times one minus hazard adjusted by replication count (N x R):
          lm_n_2 <- lm_n*(1-lm)
          
          # first and negative second derivatives of log likelihood (K x R):
          dl <- matrix(0,K[m],R) # first derivative
          nd2l <- matrix(0,K[m],R) # negative second derivative
          for (k in seq_len(K[m])) {
            subi <- subindex[[m]][[k]]
            dl[k,] <- y_sum[[m]][,k] - colSums(lm_n[subi,,drop=FALSE])
            nd2l[k,] <- colSums(lm_n_2[subi,,drop=FALSE])
          }
          
          # full conditional precision matrices (K x K x R):
          Q_fc <- array(0,c(K[m],K[m],R))
          # fill in K x K matrices separately for each risk R:
          for (r in seq_len(R)) {
            # fill in the diagonal elements:
            for (k in seq_len(K[m])) {
              Q_fc[k,k,r] <- ifelse(k==1 || k==K[m],1,2)*prec_m[[m]][r] + nd2l[k,r]
            }
            # fill in the off diagonals:
            for (k in seq_len(K[m]-1)) {
              Q_fc[k+1,k,r] <- -1*prec_m[[m]][r]
              Q_fc[k,k+1,r] <- -1*prec_m[[m]][r]
            }
          }
          
          # upper-triangular Cholesky factors of precisions (K x K x R):
          Lt_Q_fc <- array(0,c(K[m],K[m],R))
          for (r in 1:R) Lt_Q_fc[,,r] <- chol(Q_fc[,,r])
          
          # full conditional canonical b (K x R):
          cb_fc <- dl + nd2l*t(b_m[[m]])
          
          # full conditional means (K x R):
          mu_fc <- matrix(0,K[m],R)
          for (r in 1:R) mu_fc[,r] <- solve(Q_fc[,,r],cb_fc[,r])
          
          # i.i.d. standard normals for generating the proposal (K x R):
          z_star <- matrix(rnorm(K[m]*R),K[m],R)
          
          # proposed parameter update (R x K):
          b_star <- matrix(0,R,K[m])
          for (r in 1:R) b_star[r,] <- (mu_fc[,r] +
                                          backsolve(Lt_Q_fc[,,r],z_star[,r]))
          
          # step 2: calculate log proposal ratio -------------------------------
          
          # proposed step (R x K):
          prop_step <- b_star - b_m[[m]]
          
          # linear predictors under proposed values (R x N):
          eta_star <- eta + prop_step[,x[,m]]
          
          # exponentials of linear predictors under proposed values (R x N):
          e_eta_star <- exp(eta_star)
          
          # sums of exponentials of linear predictors under proposed values (N):
          cs <- .colSums(e_eta_star,R,N)
          
          # negative log partition function under proposed values (N):
          nlpartfun_star <- -log1p(cs)
          
          # discrete hazards under proposed values (N x R):
          lm_star <- t(e_eta_star)/(1+cs)
          
          # discrete hazards under proposed values adj. by replications (N X R):
          lm_n_star <- n*lm_star
          
          # hazard times one minus hazard under proposed values adj. (N x R):
          lm_n_2_star <- lm_n_star*(1-lm_star)
          
          # first and negative second derivatives of ll under proposals (K x R):
          dl_star <- matrix(0,K[m],R) # first derivative
          nd2l_star <- matrix(0,K[m],R) # negative second derivative
          for (k in seq_len(K[m])) {
            subi <- subindex[[m]][[k]]
            dl_star[k,] <- y_sum[[m]][,k] - colSums(lm_n_star[subi,,drop=FALSE])
            nd2l_star[k,] <- colSums(lm_n_2_star[subi,,drop=FALSE])
          }
          
          # full conditional precision matrices under proposals (K x K x R):
          Q_fc_star <- array(0,c(K[m],K[m],R))
          # fill in K x K matrices separately for each risk R:
          for (r in seq_len(R)) {
            # fill in the diagonal elements:
            for (k in seq_len(K[m])) {
              Q_fc_star[k,k,r] <- (ifelse(k==1 || k==K[m],1,2)*prec_m[[m]][r]
                                   + nd2l_star[k,r])
            }
            # fill in the off diagonals:
            for (k in seq_len(K[m]-1)) {
              Q_fc_star[k+1,k,r] <- -1*prec_m[[m]][r]
              Q_fc_star[k,k+1,r] <- -1*prec_m[[m]][r]
            }
          }
          
          # upper-triangular Cholesky of precisions under proposals (K x K x R):
          Lt_Q_fc_star <- array(0,c(K[m],K[m],R))
          for (r in 1:R) Lt_Q_fc_star[,,r] <- chol(Q_fc_star[,,r])
          
          # full conditional canonical b under proposals (K x R):
          cb_fc_star <- dl_star + nd2l_star*t(b_star)
          
          # full conditional means (K x R):
          mu_fc_star <- matrix(0,K[m],R)
          for (r in 1:R) mu_fc_star[,r] <- solve(Q_fc_star[,,r],cb_fc_star[,r])
        
          # calculate log proposal ratio:
          lqr <- 0
          for (r in 1:R) lqr <- (lqr - sum(log(diag(Lt_Q_fc[,,r])))
                                 + sum(log(diag(Lt_Q_fc_star[,,r])))
                                 - 0.5*crossprod(Lt_Q_fc_star[,,r] %*% (b_m[[m]][r,]
                                                -mu_fc_star[,r]))
                                 + 0.5*crossprod(z_star[,r])
          )

          # step 3: calculate log prior ratio and log likelihood ratio ---------
          
          # calculate log prior ratio:
          lpr <- 0
          for (k in seq_len(K[m]-1)) {
            lpr <- lpr + sum(dnorm(b_star[,k+1],b_star[,k],s_m[[m]],TRUE)
                             - dnorm(b_m[[m]][,k+1],b_m[[m]][,k],s_m[[m]],TRUE))
          }
          
          # calculate log likelihood ratio:
          llr <- sum(y_sum[[m]]*prop_step) + crossprod(n,nlpartfun_star-nlpartfun)
          
          # step 4: accept or reject the proposal ------------------------------
          
          # accept the proposal with the appropriate probability:
          if (log(runif(1)) < (lqr + lpr + llr)) {
            b_m[[m]] <- b_star
            eta <- eta_star
            nlpartfun <- nlpartfun_star
            if (s > B) b_m_a[[m]] <- b_m_a[[m]] + 1L
          }
                    
        }
        
        # sweep parameter vector averages onto intercepts:
        rms <- rowMeans(b_m[[m]])
        b_m[[m]] <- b_m[[m]] - rms
        b_0 <- b_0 + rms
        
        
        # update GMRF variance parameters --------------------------------------
        
        
        # calculate sums of squares of adjacent parameters for all R risks:
        SS <- apply(b_m[[m]],1, function(v) crossprod(diff(v)) )
        
        # sample the increment precisions from the full conditionals:
        prec_m[[m]] <- rgamma(R,priors[[m]]$a+0.5*(K[m]-1),priors[[m]]$b+0.5*SS)
        
        # update the corresponding std. dev. value:
        s_m[[m]] <- sqrt(1/prec_m[[m]])
                    
                
      # random effects case RWM updates ----------------------------------------
        
            
      } else if (priors[[m]]$type == "re") {
        
        for (k in seq_len(K[m])) {
          
          # current parameter values:
          b_old <- b_m[[m]][,k]
          
          # full conditional precision matrix:
          fc_prec <- prec_m[[m]] + diag(y_sum[[m]][,k])
          
          # proposed step in the random walk:
          prop_step <- 2.4/sqrt(R)*backsolve(chol(fc_prec),rnorm(R))
          
          # proposed new parameter values:
          b_star <- b_old + prop_step
          
          # log prior ratio:
          lpr <- -0.5*( (b_star %*% prec_m[[m]] %*% b_star) -
                        (b_old  %*% prec_m[[m]] %*% b_old ) )
          
          # index subset of observations with value k for variable m:
          subi <- subindex[[m]][[k]]
          
          # if there are no observations with value k for variable m, then
          # accept the proposal based on the log prior ratio only:
          if (length(subi) == 0L) {
            
            if (log(runif(1)) < lpr) {
              b_m[[m]][,k] <- b_star
              if (s > B) b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
            }

          # otherwise, calculate the log likelihood ratio then accept/reject:
          } else {
            
            # linear predictors under proposed values:
            eta_star <- eta[,subi] + prop_step
            
            # negative log partition function under proposed values:
            nlpartfun_star <- -log1p(.colSums(exp(eta_star),R,length(subi)))
            
            # log likelihood ratio:
            llr <- (crossprod(y_sum[[m]][,k],prop_step) + 
                      crossprod(n[subi],nlpartfun_star-nlpartfun[subi]))
            
            # accept the proposal with the appropriate probability:
            if (log(runif(1)) < (llr + lpr)) {
              b_m[[m]][,k] <- b_star
              eta[,subi] <- eta_star
              nlpartfun[subi] <- nlpartfun_star
              if (s > B) b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
            }
            
          }
          
        } # end k loop
        
        # update random effects covariance matrix:
        prec_m[[m]] <- rWishart(1,priors[[m]]$nu + K[m],
                                solve(priors[[m]]$s+tcrossprod(b_m[[m]])))[,,1]
          
      } # end random effects case
      
    } # end m loop
    
    
    # store scan ---------------------------------------------------------------
    
    
    # store updates only for post-burn-in iterations:
    if (s > B) {
      b_0_s[,s-B] <- b_0
      for (m in seq_len(M)) {
        if (priors[[m]]$type %in% c("cat","gmrf") ||
            (priors[[m]]$type == "re" && store_re)) {
          b_m_s[[m]][,,s-B] <- b_m[[m]]
        }
        if (priors[[m]]$type == "gmrf") s_m_s[[m]][,s-B] <- s_m[[m]]
        if (priors[[m]]$type == "re") s_m_s[[m]][,,s-B] <- solve(prec_m[[m]])
      }
    }
    
  } # end s loop
  
  # return the linear predictor params, scale params, and acceptance counters:
  list(b_0_s=b_0_s,b_m_s=b_m_s,s_m_s=s_m_s,b_m_a=b_m_a)
  
  
} # end brea_mcmc function
