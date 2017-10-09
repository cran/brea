

# This function checks that the outcome vector/matrix y and the vector n of
# counts of replicated observations are valid, and then returns integer-type
# versions of y and n in a list.


check_set_outcome <- function(y,n,N) {
  
  
  # check and transform y ------------------------------------------------------
  
  
  # make sure y is a vector or matrix of nonnegative whole numbers:
  if (!((is.vector(y) || is.matrix(y)) && all_whole(y) && all(y > -.5))) {
    stop("y must be a vector or matrix of nonnegative whole numbers")
  }
  
  # convert y from a vector to a single-column matrix if necessary:
  if (is.vector(y)) dim(y) <- c(length(y),1)
  
  # convert y from floating points to integers if necessary:
  if (!is.integer(y)) {
    y <- round(y)
    storage.mode(y) <- "integer"
  }
  
  # make sure y has N rows (same as the number of rows of x):
  if (nrow(y) != N) stop("x and y must have the same number of rows")

  # make sure there is at least one observed event for each competing risk:  
  if (any(colSums(y) < 1L)) {
    stop("there must be at least one observed event for each event type")
  }
  
  
  # check and transform n ------------------------------------------------------
  
  
  # make sure n is NULL or a vector of positive integers of the correct length:
  if (!(is.null(n) || (is.vector(n) && length(n)==N && all_whole(n)
                       && all(n > .5)))) {
    stop("n must be NULL or a vector of positive integers of length nrow(x)")
  }
  
  # if n is NULL, set it to all 1's:
  if (is.null(n)) n <- rep(1L,N)
  
  # convert n from floating points to integers if necessary:
  if (!is.integer(n)) n <- as.integer(round(n))

    
  # check compatibility of y and n and return ----------------------------------
  
  
  # make sure the number of observed events at each row is not higher than the
  # number of observations represented by that row:
  if (any(rowSums(y) > n)) stop(paste("the total number of events of all types",
                                      "cannot exceed n"))
  
  # make sure not all observations result in event occurrence:
  if (sum(y) >= sum(n)) stop("events cannot occur on all observations")

  # return the new y and n values:  
  list(y,n)
  
  
}
