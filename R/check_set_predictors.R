

# This function checks that the predictor matrix/dataframe x and the vector of
# maximum discretized covariate values K (if supplied) are valid, and then
# returns versions of x and K of type integer in a list.


check_set_predictors <- function(x,K) {
  
  
  # Case 1: x is a matrix of positive integers ---------------------------------
  
  
  if (is.matrix(x) && all_whole(x) && all(x > .5)) {
    
    # convert x from floating points to integers if necessary:
    if (!is.integer(x)) {
      x <- round(x)
      storage.mode(x) <- "integer"
    }
    
    # set K if needed, and otherwise verify K is a vector of integers:
    if (is.null(K)) {
      K <- apply(x,2,max)
    } else if (!(is.vector(K) && length(K)==ncol(x) && all_whole(K))) {
      stop("K must be NULL or a vector of whole numbers of length ncol(x)")
    }
    
    # convert K from floating points to integers if necessary:
    if (!is.integer(K)) K <- as.integer(round(K))
  
  
  # Case 2: x is a dataframe with integer or factor columns  -------------------
  
  
  } else if (is.data.frame(x)) {
    
    # first make sure K is potentially valid:
    if (!(is.null(K) || (is.vector(K) && length(K)==ncol(x) && all_whole(K)))) {
      stop("K must be NULL or a vector of whole numbers of length ncol(x)")
    }
    
    # if K is a floating point vector, convert it to integer type:
    if (!is.null(K) && !is.integer(K)) K <- as.integer(round(K))
    
    # create variables to store new x and K values:
    x.new <- matrix(0L,nrow(x),ncol(x))
    K.new <- integer(ncol(x))

    # loop over columns m of x:      
    for (m in seq_len(ncol(x))) {
      
      # use shorter variable name for x[[m]] for convenience:
      xm <- x[[m]]
      
      # case i: the predictor is a factor:
      if (is.factor(xm)) {
        
        # store integer factor codes into the new integer matrix:
        x.new[,m] <- as.integer(xm)
        
        # calculate K[m] if K is null; otherwise use supplied value:
        if (is.null(K)) {
          K.new[m] <- nlevels(xm)
        } else {
          K.new[m] <- K[m]
        }
        
      # case ii: the predictor is a vector of positive integers:
      } else if (is.vector(xm) && all_whole(xm) && all(xm > .5)) {
        
        # convert xm from floating points to integers if necessary:
        if (!is.integer(xm)) xm <- as.integer(round(xm))
        
        # store the integer-type variable xm into the new integer matrix:
        x.new[,m] <- xm
        
        # calculate K[m] if K is null; otherwise use supplied value:
        if (is.null(K)) {
          K.new[m] <- max(xm)
        } else {
          K.new[m] <- K[m]
        }
        
      # otherwise the predictor was an invalid type:
      } else stop("columns of x must be factors or positive integer vectors")
        
    } # end loop over columns m of x
    
    # replace x and K with their new values:
    x <- x.new
    K <- K.new
    
        
  # otherwise x was of neither allowed type ------------------------------------

        
  } else {
      stop(paste("x must be a matrix of positive whole numbers",
                 "or a dataframe with whole number or factor columns"))
  }
  
  
  # verify K is valid and return x and K ---------------------------------------
  
  
  # make sure there are not covariate values larger than corresponding K:      
  if (any(K < apply(x,2,max))) {
    stop("K values cannot be smaller than corresponding x variable values")
  }
  
  # make sure each covariate has (potentially) more than one value:
  if (any(K < 2L)) stop("at least one covariate has only one level")

  # return the new x and K values:
  list(x,K)
  
  
}
