

# This function checks that the prior specification is valid, and returns a
# version of the prior list with default values (if the user specified NULL for
# for the prior specification) or with additional useful forms of the prior
# parameters added.


check_set_priors <- function(priors,M,R) {
  
  
  # if priors aren't supplied, use default categorical prior for all predictors:
  if (is.null(priors)) {
    priors <- vector("list",M)
    for (m in seq_len(M)) priors[[m]] <- list("cat",4)
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
      if (length(priors[[m]]) != 2 || !is.numeric(priors[[m]][[2]]) ||
          priors[[m]][[2]] <= 0) {
        stop("positive prior full conditional std. dev. must be specified")
      }
      # set the full conditional precision value and name the parameters:
      priors[[m]][[3]] <- 1/priors[[m]][[2]]^2
      # name the parameters:
      names(priors[[m]]) <- c("type","sd","tau")
      
    # the Gausian Markov random field case:
    } else if (priors[[m]][[1]] == "gmrf") {
      if (length(priors[[m]]) != 3 || priors[[m]][[2]] <= 2 ||
          priors[[m]][[3]] <= 0) {
        stop(paste("prior degrees of freedom parameter greater than 2",
                   "and positive scale parameter must be specified"))
      }
      # set the corresponding alpha and beta for the inverse gamma distribution:
      priors[[m]][[4]] <- priors[[m]][[2]]/2
      priors[[m]][[5]] <- priors[[m]][[2]]/2*priors[[m]][[3]]
      # name the parameters:
      names(priors[[m]]) <- c("type","nu","s","a","b")
      
    # the random effects case:
    } else if (priors[[m]][[1]] == "re") {
      if (length(priors[[m]]) != 3) {
        stop("must supply prior degrees of freedom and scale matrix")
      }
      if (!is.numeric(priors[[m]][[2]]) || priors[[m]][[2]] < R + 2) {
        stop("prior degrees of freedom must be at least R + 2")
      }
      if (!(is.matrix(priors[[m]][[3]]) && all(dim(priors[[m]][[3]]) == R) &&
            isSymmetric(priors[[m]][[3]]) &&
            all(eigen(priors[[m]][[3]],TRUE)$values > 0))) {
        stop("must supply positive-definite R x R prior scale matrix")
      }
      # name the parameters:
      names(priors[[m]]) <- c("type","nu","s")
      
    # the zero (unused) case:
    } else if (priors[[m]][[1]] == "zero") {
      names(priors[[m]])[1] <- "type"
    
    # otherwise, the prior was none of the allowed types:
    } else stop("prior type must be cat, gmrf, re, or zero")
    
  } # end m loop
  
  # return the updated prior specs:
  priors
  
  
}
