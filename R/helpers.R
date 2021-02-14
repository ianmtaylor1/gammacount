# This file is for functions which:
#  1. Perform relatively basic tasks
#  2. Are used by potentially very many functions in other files
#  3. Don't have a logical home elsewhere.
# These functions can be thought of as functions that could be useful additions
# to R generally.

# Helper function to compute log(exp(a) - exp(b))
# a and b should be of equal length, but no checking is done
# There are seventeen total cases for the input values of this function.
# a and b can both be either: NaN, -Inf, Inf, or finite. If both are finite,
# Then either a >= b or b < a.
# In all cases, NaN should result if appropriate:
#   1. Either a or b is NaN
#   2. b > a (b either finite or infinite)
#   3. Both a and b are +Inf
# Every other case should result in a non-NaN value
logdiffexp <- function(a, b) {

  diffs <- b - a

  # Special cases
  bothneginf <- is.infinite(a) & is.infinite(b) & (a == -Inf) & (b == -Inf)
  smalldiff <- is.finite(a) & is.finite(b) & (abs(diffs) < log(2))
  defaultcase <- (!bothneginf) & (!smalldiff)

  # Split computation according to
  # https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
  val <- rep(NaN, length(diffs))
  # The general case
  val[defaultcase] <- a[defaultcase] + log1p(-exp(diffs[defaultcase]))
  # The small difference case
  val[smalldiff] <- a[smalldiff] + log(-expm1(diffs[smalldiff]))
  # The negative infinity case
  val[bothneginf] <- -Inf

  # Alert of any new NaNs
  # Everything except the positive infinity case should be caught by log
  # or log1p
  if (any(is.infinite(a) & is.infinite(b) & (a > 0) & (b > 0))) {
    warning("NaNs produced")
  }

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
