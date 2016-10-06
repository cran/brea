



brea_mcmc <- function(x, y, S = 1000L, priors = NULL, n = NULL, K = NULL,
                      store_re = FALSE) {

  
  # verify and set x and K -----------------------------------------------------

  
  # if given a positive integer matrix x, set the number of predictor levels K
  # if not supplied; otherwise, check that K is a positive integer vector of the
  # correct length:
  if (is.matrix(x) && is.integer(x) && all(x > 0)) {
    if (is.null(K)) {
      K <- apply(x,2,max)
    } else if (!is.vector(K) || length(K) != ncol(x)) {
      stop("K must be a vector with length equal to the number of columns of x")
    } else if (!is.integer(K) || any(K < 1)) {
      stop("K must consist of positive integers")
    }
  
  # if given a dataframe x with all columns factors, set the number of predictor
  # levels K and then convert x to an integer matrix:
  } else if (is.data.frame(x) && all(sapply(x,is.factor))) {
    K <- sapply(x,nlevels)
    x <- sapply(x,as.integer)
    
  # otherwise the supplied x is neither allowed type:
  } else stop("x must be a positive integer matrix or factor dataframe")
  
  # make sure each covariate assumes more than one value:
  if (any(K < 2)) stop("at least one covariate has only one level")
  # note that this does not prevent the user from supplying a dataset in which
  # some covariates only assume one value
  
  # set the number of predictors:
  M <- ncol(x)
  
  # set the number of observations:
  N <- nrow(x)

  
  # verify y is valid ----------------------------------------------------------

  
  # if y is a vector, make it a column vector:
  if (is.vector(y)) y <- cbind(y)
  
  # make sure the number of rows of y and x are equal:
  if (nrow(y) != N) stop("x and y must have the same number of rows")
  
  # make sure y consists of nonnegative integers:
  if (!is.integer(y) || any(y < 0)) {
    stop("y must be an integer matrix with nonnegative entries")
  }
  
  # make sure there is at least one observed event for each competing risk: 
  if (any(colSums(y) < 1)) {
    stop("there must be at least one observed event for each event type")
  }
  
  # set the number of event types (i.e. competing risks):
  R <- ncol(y)

    
  # verify n is valid ----------------------------------------------------------

  
  # set the number of observations each row represents if not supplied:
  if (is.null(n)) {
    n <- rep(1L,N)
  # otherwise check that the number of observations is valid:
  } else if (!is.integer(n) || any(n < 1)) {
    stop("n must be an integer vector with positive entries")
  # check that number of observed events does not exceed number of observations:
  } else if (any(rowSums(y) > n)) {
    stop("the total number of events of all types cannot exceed n")
  }

    
  # verify prior type & specification are valid --------------------------------


  # if priors aren't supplied, use default categorical prior for all predictors:
  if (is.null(priors)) {
    priors <- vector("list",M)
    for (m in seq_len(M)) priors[[m]] <- list("cat",2)
  }

  # check that priors is a list of the appropriate length:
  if (!is.list(priors) || length(priors) != M) {
    stop("priors must be a list with one element for each column of x")
  }
  
  # check that each element is a valid prior specification, and name the
  # parameters of each prior for convenience:
  for (m in seq_len(M)) {
    
    # make sure prior spec is a list of length at least 1:
    if (!is.list(priors[[m]]) || length(priors[[m]]) == 0) {
      stop("each prior specification must be a list with at least one element")
      
    # the categorical case:
    } else if (priors[[m]][[1]] == "cat") {
      if (length(priors[[m]]) < 2 || !is.numeric(priors[[m]][[2]]) ||
          priors[[m]][[2]] <= 0) {
        stop("positive prior full conditional std. dev. must be specified")
      }
      # set the full conditional precision value and name the parameters:
      priors[[m]][[3]] <- 1/priors[[m]][[2]]^2
      names(priors[[m]])[2:3] <- c("sd","tau")
    
    # the Gausian Markov random field case:
    } else if (priors[[m]][[1]] == "gmrf") {
      if (length(priors[[m]]) < 3 || priors[[m]][[2]] <= 0 ||
          priors[[m]][[3]] <= 0) {
        stop("positive degrees of freedom and scale must be specified")
      }
      # set the corresponding alpha and beta for the inverse gamme distribution:
      priors[[m]][[4]] <- priors[[m]][[2]]/2
      priors[[m]][[5]] <- priors[[m]][[2]]/2*priors[[m]][[3]]
      # name the parameters:
      names(priors[[m]])[2:5] <- c("nu","s","a","b")

    # the random effects case:
    } else if (priors[[m]][[1]] == "re") {
      if (length(priors[[m]]) < 3) {
        stop("must supply prior degrees of freedom and scale matrix")
      }
      if (priors[[m]][[2]] < R + 2) {
        stop("prior degrees of freedom must be at least R + 2")
      }
      if (!is.matrix(priors[[m]][[3]]) || !all(dim(priors[[m]][[3]]) == R)) {
        stop("must supply R x R prior scale matrix")
      }
      # name the parameters:
      names(priors[[m]])[2:3] <- c("nu","s")

    # if none of the the above or not the unused (zero) case, fail:
    } else if (priors[[m]][[1]] != "zero") {
      stop("prior type must be cat, gmrf, re, or zero")
    }
  
  } # end m loop
  

  # initialize parameters, parameter storage, and acceptance counters ----------

  
  # check that the number S of requested scans is a positive integer:
  if (!is.numeric(S) || (S %% 1 != 0) || (S < 1)) {
    stop("S must be a nonnegative integer")
  }
  
  # initialize intercepts to log-odds of event occurrence for each risk:
  logit <- function(p) log(p/(1-p))
  b_0 <- logit(colSums(y)/sum(n))
  
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
    if (priors[[m]][[1]] %in% c("cat","gmrf") ||
        (priors[[m]][[1]] == "re" && store_re)) {
      b_m_s[[m]] <- array(0,c(R,K[m],S))
    }
    # acceptance counters:
    b_m_a[[m]] <- rep(0L,K[m])
  }
  
  # current values, storage and acceptance counters for variance/precision
  # hyperparameters for the GMRF and random effects cases:
  
  # std. dev. of the random walk increment or covariance of the random effects:
  s_m <- vector("list",M)
  s_m_s <- vector("list",M)
  # precision of the random walk increment or precision of the random effects:
  prec_m <- vector("list",M)

  for (m in seq_len(M)) {
    
    # GMRF case:
    if (priors[[m]][[1]] == "gmrf") {
      # initialize increment std. dev. to 1:
      s_m[[m]] <- rep(1,R)
      s_m_s[[m]] <- matrix(0,R,S)
      prec_m[[m]] <- 1/s_m[[m]]^2
      
    # random effects case:
    } else if (priors[[m]][[1]] == "re") {
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
  for (m in seq_len(m)) eta <- eta + b_m[[m]][,x[,m]]
  
  # current values of the negative log partition functions:
  nlpartfun <- -log1p(colSums(exp(eta)))
  
  
  # create index subsets and y-sums for single-category RWM steps --------------
  
  
  # list whose m^th element is the list whose k^th element is the vector of
  # indices of observations with value k for predictor m:
  subindex <- vector("list",M)
  
  # list whose m^th element is the R x K[m] matrix whose (r,k) entry is the
  # count of events of type r occurring at observatsion with value k for
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

  
  # loop over scans s:  
  for (s in seq_len(S)) {
  
    # loop over predictors m:  
    for (m in seq_len(M)) {
      
      
      # categorical case RWM updates -------------------------------------------
      
      
      if (priors[[m]][[1]] == "cat") {
        
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
          
          # index subset with value k for variable m:
          subi <- subindex[[m]][[k]]
          
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
            b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
          }
          
        } # end k loop
        
        # sweep parameter vector averages onto intercepts:
        rms <- rowMeans(b_m[[m]])
        b_m[[m]] <- b_m[[m]] - rms
        b_0 <- b_0 + rms
        
      
      # GMRF case RWM updates --------------------------------------------------
      
        
      } else if (priors[[m]][[1]] == "gmrf") {
        
        for (k in seq_len(K[m])) {
          
          # full conditional precision:
          fc_prec <- ifelse(k==1 || k==K[m],1,2)*prec_m[[m]] + y_sum[[m]][,k]
          
          # proposed step in the random walk:
          prop_step <- rnorm(R,0,2.4/sqrt(R*fc_prec))
          
          # proposed new parameter values:
          b_star <- b_m[[m]][,k] + prop_step
          
          # log prior ratio:
          if (k==1) {
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
          
          # index subset with value k for variable m:
          subi <- subindex[[m]][[k]]
          
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
            b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
          }
          
        } # end k loop
        
        # sweep parameter vector averages onto intercepts:
        rms <- rowMeans(b_m[[m]])
        b_m[[m]] <- b_m[[m]] - rms
        b_0 <- b_0 + rms
        
        
        # update GMRF hyperparameters ------------------------------------------
        
        
        # calculate sums of squares of adjacent parameters for all R risks:
        SS <- apply(b_m[[m]],1, function(v) crossprod(diff(v)) )
        
        # sample the increment precisions from the full conditionals:
        prec_m[[m]] <- rgamma(R,priors[[m]]$a+.5*(K[m]-1),priors[[m]]$b+.5*SS)
        
        # update the corresponding std. dev. value:
        s_m[[m]] <- sqrt(1/prec_m[[m]])
                    
                
      # random effects case RWM updates ----------------------------------------
        
            
      } else if (priors[[m]][[1]] == "re") {
        
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
          
          # index subset with value k for variable m:
          subi <- subindex[[m]][[k]]
          
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
            b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
          }
          
        } # end k loop
        
        # update random effects covariance matrix:
        prec_m[[m]] <- rWishart(1,priors[[m]]$nu + K[m],
                                solve(priors[[m]]$s+tcrossprod(b_m[[m]])))[,,1]
          
      } # end random effects case
      
    } # end m loop
    
    
    # store scan ---------------------------------------------------------------
    
    
    b_0_s[,s] <- b_0
    for (m in seq_len(M)) {
      if (priors[[m]][[1]] %in% c("cat","gmrf") ||
          (priors[[m]][[1]] == "re" && store_re)) {
        b_m_s[[m]][,,s] <- b_m[[m]]
      }
      if (priors[[m]][[1]] == "gmrf") s_m_s[[m]][,s] <- s_m[[m]]
      if (priors[[m]][[1]] == "re") s_m_s[[m]][,,s] <- solve(prec_m[[m]])
    }
    
  } # end s loop
  
  # return the linear predictor params, scale params, and acceptance counters:
  list(b_0_s=b_0_s,b_m_s=b_m_s,s_m_s=s_m_s,b_m_a=b_m_a)
  
} # end brea_mcmc function

