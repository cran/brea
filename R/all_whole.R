

# This function tests whether the passed object x consists of a vector/array of
# whole numbers (either integer or floating point type). It is based upon the
# is.wholenumber function found in the is.integer help file.


all_whole <- function(x, tol = sqrt(.Machine$double.eps)) {
  is.integer(x) || (is.numeric(x) && all(abs(x - round(x)) < tol))
}
