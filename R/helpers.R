# This file is for functions which:
#  1. Perform relatively basic tasks
#  2. Are used by potentially very many functions in other files
#  3. Don't have a logical home elsewhere.
# These functions can be thought of as functions that could be useful additions
# to R generally.

# Helper function to compute log(exp(a) - exp(b))
# a and b should be of equal length, but no checking is done
logdiffexp <- function(a, b) {
  if (any(a < b)) {
    warning("a must be larger than b. NaNs produced")
  }
  # Split computation according to
  # https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
  diffs <- b - a
  diffs[(b == -Inf) & (a == -Inf)] <- 0 # Fix case when a == b == -Inf
  smalldiff <- (abs(diffs) < log(2))
  val <- rep(0, length(diffs))
  # For large differences, use log1p(-exp())
  val[!smalldiff] <- log1p(-exp(diffs[!smalldiff]))
  # For small differences, use log(-expm1())
  val[smalldiff] <- log(-expm1(diffs[smalldiff]))
  # Add a to all values
  val <- a + val
  return(val)
}

# Helper function to compute log(exp(a) + exp(b)). Useful for adding
# probabilities on the log scale
logsumexp <- function(a, b) {
  negdiff <- -abs(a - b)
  negdiff[(a == -Inf) & (b == -Inf)] <- 0 # Fix cases where a == b == -Inf
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
