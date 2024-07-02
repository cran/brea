
brea_mcmc_ms <- function(x, y, priors = NULL, S = 1000, B = 100, n = NULL,
                      K = NULL, store_re = FALSE) {

  
  # verify and clean user input ------------------------------------------------
  
  # check that x and y are lists of the same nonzero length:
  if (!(is.list(x) && is.list(y) && length(x)==length(y) && length(x) > 0L))
    stop("x and y must be lists of the same nonzero length")
  
  # set the number of states:
  ST <- length(x)
  
  # if priors is NULL, make it a list of NULL's;
  # otherwise check that it's a list of length ST:
  if (is.null(priors)) {
    priors <- vector("list",ST)
  } else if (!(is.list(priors) && length(priors)==ST)) {
    stop("priors must be NULL or a list of the same length as x and y")
  }
  
  # if n is NULL, make it a list of NULL's;
  # otherwise check that it's a list of length ST:
  if (is.null(n)) {
    n <- vector("list",ST)
  } else if (!(is.list(n) && length(n)==ST)) {
    stop("n must be NULL or a list of the same length as x and y")
  }
  
  # if K is NULL, make it a list of NULL's;
  # otherwise check that it's a list of length ST:
  if (is.null(K)) {
    K <- vector("list",ST)
  } else if (!(is.list(K) && length(K)==ST)) {
    stop("K must be NULL or a list of the same length as x and y")
  }
  
  # verify and set x[[st]] and K[[st]] for each state st:
  for (st in 1:ST) {
    xK <- check_set_predictors(x[[st]],K[[st]])
    x[[st]] <- xK[[1]]
    K[[st]] <- xK[[2]]
  }
  
  # set the number of predictors for each state:
  M <- sapply(x,ncol)
  
  # set the number of observations for each state:
  N <- sapply(x,nrow)
  
  # verify and set y[[st]] and n[[st]] for each state st:
  for (st in 1:ST) {
    yn <- check_set_outcome(y[[st]],n[[st]],N[st])
    y[[st]] <- yn[[1]]
    n[[st]] <- yn[[2]]
  }
  
  # set the number of event types (i.e. competing risks):
  R <- sapply(y,ncol)
  
  # verify and set priors[[st]] for each state st:
  for (st in 1:ST) {
    priors[[st]] <- check_set_priors(priors[[st]],M[st],sum(R))
  }
  
  
  # verify that "re" prior type indices m agree for all states -----------------
  
  
  # make sure there are no "re" priors for m > min(M):
  for (st in 1:ST) {
    if (M[st] > min(M)) for (m in (min(M)+1):M[st]) {
      if (priors[[st]][[m]]$type == "re") {
        stop("columns of x assigned re prior types must agree across states")
      }
    }
  }
  
  # for m <= min(M), determine which m have "re" type for each state:
  re_indices <- matrix(0L,ST,min(M))
  for (st in 1:ST) for (m in 1:min(M)) {
    if (priors[[st]][[m]]$type == "re") re_indices[st,m] <- 1L
  }
  
  # check if states disagree on re prior type:
  if (!all(colSums(re_indices) %in% c(0,ST))) {
    stop("columns of x assigned re prior types must agree across states")
  }

  # indices of components of the joint sum(R)-length random effects vector
  # belonging to each state:
  ind_st <- vector("list",ST)
  ni <- 1 # next available index
  for (st in 1:ST) {
    ind_st[[st]] <- ni:(ni+R[st]-1)
    ni <- ni + R[st]
  }
 
  # check whether K values agree across states for all re indices:
  for (m in which(re_indices[1,]==1L)) {
    K_st <- integer(ST)
    for (st in 1:ST) K_st[st] <- K[[st]][m]
    # warn if not all K values are the same, and then set them all to max(K_st):
    if (!all(diff(K_st)==0)) {
      warning(paste("K values do not agree across states for re index",m))
      for (st in 1:ST) K[[st]][m] <- max(K_st)
    }
  }
  
  
  # verify correctness of run length -------------------------------------------
  
    
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

  
  # linear predictor parameters:
  b_0 <- vector("list",ST)
  b_0_s <- vector("list",ST)
  b_m <- vector("list",ST)
  b_m_s <- vector("list",ST)
  b_m_a <- vector("list",ST)
  
  for (st in 1:ST) {
    
    # initialize current values of intercepts to MLE's with covariates omitted:
    b_0[[st]] <- log(colSums(y[[st]])/(sum(n[[st]])-sum(y[[st]])))
    
    # intercept storage:
    b_0_s[[st]] <- matrix(0,R[st],S)
    
    # current values, storage, and acceptance counters for linear pred. params.:
    b_m[[st]] <- vector("list",M[st])
    b_m_s[[st]] <- vector("list",M[st])
    b_m_a[[st]] <- vector("list",M[st])
    
    for (m in seq_len(M[st])) {
      # current parameter values:
      b_m[[st]][[m]] <- matrix(0,R[st],K[[st]][m])
      # stored parameter values (saving random effects only if requested):
      if (priors[[st]][[m]]$type %in% c("cat","gmrf") ||
          (priors[[st]][[m]]$type == "re" && store_re)) {
        b_m_s[[st]][[m]] <- array(0,c(R[st],K[[st]][m],S))
      }
      # acceptance counters:
      b_m_a[[st]][[m]] <- rep(0L,K[[st]][m])
    }
    
  }
  
  # variance/precision parameters:
  s_m <- vector("list",ST)
  s_m_s <- vector("list",ST)
  prec_m <- vector("list",ST)
  
  for (st in 1:ST) {
    
    # std. dev. of the random walk increment or covariance of the random effects:
    s_m[[st]] <- vector("list",M[st])
    s_m_s[[st]] <- vector("list",M[st])
    
    # precision of the random walk increment or precision of the random effects:
    prec_m[[st]] <- vector("list",M[st])
    
    for (m in seq_len(M[st])) {
      
      # GMRF case:
      if (priors[[st]][[m]]$type == "gmrf") {
        # initialize increment std. dev. to 1:
        s_m[[st]][[m]] <- rep(1,R[st])
        s_m_s[[st]][[m]] <- matrix(0,R[st],S)
        prec_m[[st]][[m]] <- 1/s_m[[st]][[m]]^2
        
      # random effects case (only storing values for first state):
      } else if (priors[[st]][[m]]$type == "re" && st==1) {
        # initialize random effects covariance to identity matrix:
        s_m[[st]][[m]] <- diag(sum(R))
        s_m_s[[st]][[m]] <- array(0,c(sum(R),sum(R),S))
        prec_m[[st]][[m]] <- solve(s_m[[st]][[m]])
      }
      
    }
    
  }
  
  # initialize linear predictors and negative log partition function:
  
  # current values of the linear predictors:
  eta <- vector("list",ST)
  # negative log partition function values:
  nlpartfun <- vector("list",ST)
  
  for (st in 1:ST) {
    eta[[st]] <- matrix(0,R[st],N[st])
    eta[[st]] <- eta[[st]] + b_0[[st]]
    for (m in seq_len(M[st])) {
      eta[[st]] <- eta[[st]] + b_m[[st]][[m]][,x[[st]][,m]]
    }
    nlpartfun[[st]] <- -log1p(colSums(exp(eta[[st]])))
  }
  
  
  # create index subsets and y-sums for single-category RWM steps --------------
  
  
  subindex <- vector("list",ST)
  y_sum <- vector("list",ST)
  
  for (st in 1:ST) {
    # list whose m^th element is the list whose k^th element is the vector of
    # indices of observations with value k for predictor m:
    subindex[[st]] <- vector("list",M[st])
    
    # list whose m^th element is the R x K[m] matrix whose (r,k) entry is the
    # count of events of type r occurring at observations with value k for
    # predictor m:
    y_sum[[st]] <- vector("list",M[st])
    
    # loop over predictors m:
    for (m in seq_len(M[st])) {
      subindex[[st]][[m]] <- vector("list",K[[st]][m])
      y_sum[[st]][[m]] <- matrix(0,R[st],K[[st]][m])
      for (k in seq_len(K[[st]][m])) {
        subindex[[st]][[m]][[k]] <- which(x[[st]][,m]==k)
        y_sum[[st]][[m]][,k] <- colSums(y[[st]][subindex[[st]][[m]][[k]],
                                                ,drop=FALSE])
      }
    }
  }
  

  # run MCMC -------------------------------------------------------------------

  
  print("Setup complete...starting MCMC.")
  
  # loop over scans s from both burn-in and stored iterations:  
  for (s in seq_len(S+B)) {
    
    # loop over states st:
    for (st in seq_len(ST)) {
  
    # loop over predictors m:  
    for (m in seq_len(M[st])) {
      
      
      # categorical case RWM updates -------------------------------------------
      
      
      if (priors[[st]][[m]]$type == "cat") {
        
        for (k in seq_len(K[[st]][m])) {
          
          # proposed step in the random walk:
          prop_step <- rnorm(R[st],0,2.4/sqrt(R[st]*(
            priors[[st]][[m]]$tau + y_sum[[st]][[m]][,k])))
          
          # proposed new parameter values:
          b_star <- b_m[[st]][[m]][,k] + prop_step
          
          # conditional prior means:
          cm <- rowMeans(b_m[[st]][[m]][,-k,drop=FALSE])
          
          # calculate log prior ratio:
          lpr <- sum(dnorm(b_star,cm,priors[[st]][[m]]$sd,TRUE)
                     - dnorm(b_m[[st]][[m]][,k],cm,priors[[st]][[m]]$sd,TRUE))
          
          # index subset of observations with value k for variable m:
          subi <- subindex[[st]][[m]][[k]]
          
          # if there are no observations with value k for variable m, then
          # accept the proposal based on the log prior ratio only:
          if (length(subi) == 0L) {
            
            if (log(runif(1)) < lpr) {
              b_m[[st]][[m]][,k] <- b_star
              if (s > B) b_m_a[[st]][[m]][k] <- b_m_a[[st]][[m]][k] + 1L
            }

          # otherwise, calculate the log likelihood ratio then accept/reject:
          } else {
            
            # linear predictors under proposed values:
            eta_star <- eta[[st]][,subi] + prop_step
            
            # negative log partition function under proposed values:
            nlpartfun_star <- -log1p(.colSums(exp(eta_star),R[st],length(subi)))
            
            # log likelihood ratio:
            llr <- (crossprod(y_sum[[st]][[m]][,k],prop_step) + crossprod(
                    n[[st]][subi],nlpartfun_star-nlpartfun[[st]][subi]))
            
            # accept the proposal with the appropriate probability:
            if (log(runif(1)) < (llr + lpr)) {
              b_m[[st]][[m]][,k] <- b_star
              eta[[st]][,subi] <- eta_star
              nlpartfun[[st]][subi] <- nlpartfun_star
              if (s > B) b_m_a[[st]][[m]][k] <- b_m_a[[st]][[m]][k] + 1L
            }
            
          }
          
        } # end k loop
        
        # sweep parameter vector averages onto intercepts:
        rms <- rowMeans(b_m[[st]][[m]])
        b_m[[st]][[m]] <- b_m[[st]][[m]] - rms
        b_0[[st]] <- b_0[[st]] + rms
        
      
      # GMRF case RWM updates --------------------------------------------------
      
        
      } else if (priors[[st]][[m]]$type == "gmrf") {
        
        for (k in seq_len(K[[st]][m])) {
          
          # full conditional precision:
          fc_prec <- (ifelse(k==1L || k==K[[st]][m],1,2)*prec_m[[st]][[m]]
                      + y_sum[[st]][[m]][,k])
          
          # proposed step in the random walk:
          prop_step <- rnorm(R[st],0,2.4/sqrt(R[st]*fc_prec))
          
          # proposed new parameter values:
          b_star <- b_m[[st]][[m]][,k] + prop_step
          
          # set rw s.d. for convenience:
          sd <- s_m[[st]][[m]]
          
          # log prior ratio:
          if (k==1L) {
            lpr <- sum(dnorm(b_star,b_m[[st]][[m]][,2L],sd,TRUE) -
                       dnorm(b_m[[st]][[m]][,k],b_m[[st]][[m]][,2L],sd,TRUE))
          } else if (k==K[[st]][m]) {
            lpr <- sum(dnorm(b_star,b_m[[st]][[m]][,k-1L],sd,TRUE) -
                       dnorm(b_m[[st]][[m]][,k],b_m[[st]][[m]][,k-1L],sd,TRUE))
          } else {
            b_near <- (b_m[[st]][[m]][,k-1L] + b_m[[st]][[m]][,k+1L])/2
            lpr <- sum(dnorm(b_star,b_near,sd/sqrt(2),TRUE) -
                       dnorm(b_m[[st]][[m]][,k],b_near,sd/sqrt(2),TRUE))
          }
          
          # index subset of observations with value k for variable m:
          subi <- subindex[[st]][[m]][[k]]
          
          # if there are no observations with value k for variable m, then
          # accept the proposal based on the log prior ratio only:
          if (length(subi) == 0L) {
            
            if (log(runif(1)) < lpr) {
              b_m[[st]][[m]][,k] <- b_star
              if (s > B) b_m_a[[st]][[m]][k] <- b_m_a[[st]][[m]][k] + 1L
            }

          # otherwise, calculate the log likelihood ratio then accept/reject:
          } else {
            
            # linear predictors under proposed values:
            eta_star <- eta[[st]][,subi] + prop_step
            
            # negative log partition function under proposed values:
            nlpartfun_star <- -log1p(.colSums(exp(eta_star),R[st],length(subi)))
            
            # log likelihood ratio:
            llr <- (crossprod(y_sum[[st]][[m]][,k],prop_step) + crossprod(
                    n[[st]][subi],nlpartfun_star-nlpartfun[[st]][subi]))
            
            # accept the proposal with the appropriate probability:
            if (log(runif(1)) < (llr + lpr)) {
              b_m[[st]][[m]][,k] <- b_star
              eta[[st]][,subi] <- eta_star
              nlpartfun[[st]][subi] <- nlpartfun_star
              if (s > B) b_m_a[[st]][[m]][k] <- b_m_a[[st]][[m]][k] + 1L
            }
            
          }
          
        } # end k loop
        
        # sweep parameter vector averages onto intercepts:
        rms <- rowMeans(b_m[[st]][[m]])
        b_m[[st]][[m]] <- b_m[[st]][[m]] - rms
        b_0[[st]] <- b_0[[st]] + rms
        
        
        # update GMRF variance parameters --------------------------------------
        
        
        # calculate sums of squares of adjacent parameters for all R risks:
        SS <- apply(b_m[[st]][[m]],1, function(v) crossprod(diff(v)) )
        
        # sample the increment precisions from the full conditionals:
        prec_m[[st]][[m]] <- rgamma(R[st],priors[[st]][[m]]$a+.5*(K[[st]][m]-1),
                                    priors[[st]][[m]]$b+.5*SS)
        
        # update the corresponding std. dev. value:
        s_m[[st]][[m]] <- sqrt(1/prec_m[[st]][[m]])
        
      } # end if/else for prior types
        
    } # end m loop
      
    } # end st loop        
                    
                
    # random effects case RWM updates ----------------------------------------
       
    
    # find predictors representing random effects:
    for (m in 1:R[1]) if (priors[[1]][[m]]$type == "re") {
      
      # current random effects vector for all states and risks:
      b <- numeric(sum(R))
      
      # y sums for all states and risks:
      ys <- numeric(sum(R))
      
      # observation subsets with value k for each state:
      subi <- vector("list",ST)
      
      # proposed linear predictor values for each state:
      eta_star <- vector("list",ST)
      
      # negative log partition function for each state:
      nlpartfun_star <- vector("list",ST)
      
      # update all sum(R) random effects jointly for each category k separately:
      for (k in seq_len(K[[1]][m])) {
        
        # initialize current parameter values and y sums for category k:
        for (st in 1:ST) {
          b[ind_st[[st]]] <- b_m[[st]][[m]][,k]
          ys[ind_st[[st]]] <- y_sum[[st]][[m]][,k]
        }
        
        # full conditional precision matrix:
        fc_prec <- prec_m[[1]][[m]] + diag(ys)
        
        # proposed step in the random walk:
        prop_step <- 2.4/sqrt(sum(R))*backsolve(chol(fc_prec),rnorm(sum(R)))
        
        # proposed new parameter values:
        b_star <- b + prop_step
        
        # log prior ratio:
        lpr <- -0.5*( (b_star %*% prec_m[[1]][[m]] %*% b_star) -
                        (b %*% prec_m[[1]][[m]] %*% b) )
        
        # log likelihood ratio:
        llr <- crossprod(ys,prop_step)
        for (st in 1:ST){
          subi[[st]] <- subindex[[st]][[m]][[k]]
          if (length(subi[[st]]) > 0) {
            eta_star[[st]] <- eta[[st]][,subi[[st]]] + prop_step[ind_st[[st]]]
            nlpartfun_star[[st]] <- -log1p(.colSums(exp(eta_star[[st]]),R[st],
                                                    length(subi[[st]])))
            llr <- llr + crossprod(n[[st]][subi[[st]]],
                                   nlpartfun_star[[st]]-nlpartfun[[st]][subi[[st]]])
          }
        }

        # accept the proposal with the appropriate probability:
        if (log(runif(1)) < (llr+lpr)) for (st in 1:ST) {
          b_m[[st]][[m]][,k] <- b_star[ind_st[[st]]]
          eta[[st]][,subi[[st]]] <- eta_star[[st]]
          nlpartfun[[st]][subi[[st]]] <- nlpartfun_star[[st]]
          if (s > B) b_m_a[[st]][[m]][k] <- b_m_a[[st]][[m]][k] + 1
        }
        
      } # end k loop
      
      
      # update random effects covariance matrix --------------------------------
      
      
      # store all random effects in a single matrix:
      re <- matrix(0,sum(R),K[[1]][m])
      for (st in 1:ST) re[ind_st[[st]],] <- b_m[[st]][[m]]
      
      # sample the Wishart full conditional for the precision matrix:
      prec_m[[1]][[m]] <- rWishart(1,priors[[1]][[m]]$nu + K[[1]][m],
                                solve(priors[[1]][[m]]$s + tcrossprod(re)))[,,1]
      
    } # end random effects m loop
    
    
    # store scan ---------------------------------------------------------------
    
 
    # store updates only for post-burn-in iterations:
    if (s > B) {
      for (st in 1:ST) {
        b_0_s[[st]][,s-B] <- b_0[[st]]
        for (m in seq_len(M[st])) {
          if (priors[[st]][[m]]$type %in% c("cat","gmrf") ||
              (priors[[st]][[m]]$type == "re" && store_re)) {
            b_m_s[[st]][[m]][,,s-B] <- b_m[[st]][[m]]
          }
          if (priors[[st]][[m]]$type == "gmrf") {
            s_m_s[[st]][[m]][,s-B] <- s_m[[st]][[m]]
          }
          if (priors[[st]][[m]]$type == "re" && st==1) {
            s_m_s[[st]][[m]][,,s-B] <- solve(prec_m[[st]][[m]])
          }
        }
      }
    }
    
  } # end s loop
  
  # return the linear predictor params, scale params, and acceptance counters:
  list(b_0_s=b_0_s,b_m_s=b_m_s,s_m_s=s_m_s,b_m_a=b_m_a)
  
  
} # end brea_mcmc function
