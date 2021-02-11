# This file is for functions which:
#  1. Perform relatively basic tasks
#  2. Are used by potentially very many functions in other files
#  3. Don't have a logical home elsewhere.
# These functions can be thought of as functions that could be useful additions
# to R generally.

# Helper function to compute log(exp(a) - exp(b))
# a and b should be of equal length, but no checking is done
logdiffexp <- function(a, b) {

  diffs <- b - a

  # Define cases
  eithernan <- is.nan(a) | is.nan(b)
  bothinf <- is.infinite(a) & is.infinite(b) & (a == b)
  blarger <- (!eithernan) & (a < b)
  smalldiff <- (!eithernan) & (!bothinf) & (!blarger) & (abs(diffs) < log(2))
  largediff <- (!eithernan) & (!bothinf) & (!blarger) & (abs(diffs) >= log(2))

  if (any(blarger)) {
    warning("a must be larger than b. NaNs produced")
  }

  # Split computation according to
  # https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
  val <- rep(0, length(diffs))
  # For large differences, use log1p(-exp())
  val[largediff] <- log1p(-exp(diffs[largediff]))
  # For small differences, use log(-expm1())
  val[smalldiff] <- log(-expm1(diffs[smalldiff]))
  # For either value being nan, the result is nan
  val[eithernan] <- NaN
  # If both are infinite and equal, the result is the same value as both
  val[bothinf] <- 0
  # In the cases where b > a, the result is NaN
  val[blarger] <- NaN

  # Add a to all values
  val <- a + val

  return(val)
}

# Helper function to compute log(exp(a) + exp(b)). Useful for adding
# probabilities on the log scale
logsumexp <- function(a, b) {
  # If either a or b is NaN, that will carry through without issue
  bothinf <- is.infinite(a) & is.infinite(b) & (a == b)
  negdiff <- -abs(a - b)
  negdiff[bothinf] <- 0 # Fix cases where a == b == -Inf or Inf
  pmax(a, b) + log1p(exp(negdiff))
}

# Recycle a vector v until it is of length exactly n. Elements of v will repeat
# until there are n of them, then cut off the rest. If length(v) > n, the result
# will be shortened
# Ex: recycle(c(1,2,3), 7) produces c(1,2,3,1,2,3,1)
# Ex: recycle(c(1,2,3), 2) produces c(1,2)
recycle <- function(v, n) {
  rep(v, ceiling(n / length(v)))[1:n]
}
