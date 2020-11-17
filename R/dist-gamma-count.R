# This file contains functions for the Gamma count distribution
# Like poisson, but arrivals are gamma(alpha,alpha). (expected value = 1)
# Length of count interval is lambda, like poisson.
# poisson is a special case with alpha=1.


# Cumulative distribution function for gamma-count distribution
pgc <- function(x, lambda, alpha=1, lower.tail=TRUE, log.p=FALSE) {
  # Use pgamma to calculate the probabilities (following Zeviani et al. (2013))
  # A few useful tricks:
  #  1. Fill invalid parameter values with NaN so that the NaN propogates to the
  #     result.
  #  2. pgamma with a zero shape parameter is kind of a point mass at zero, so
  #     x == -1 behaves nicely for all x < 0. Use pmax to provide lower bound.
  pgamma(ifelse(lambda > 0, lambda, NaN),
         (pmax(floor(x), -1) + 1) * ifelse(alpha > 0, alpha, NaN),
         rate = ifelse(alpha > 0, alpha, NaN),
         lower.tail = !lower.tail,
         log.p = log.p)
}

# Probability mass funtion for gamma-count distribution
dgc <- function(x, lambda, alpha=1, log=FALSE) {
  # Define the pmf in terms of the cdf. This avoids double-checking parameter
  # values.
  lp <- logdiffexp(
    pgc(x, lambda, alpha, log.p=TRUE),
    pgc(x-1, lambda, alpha, log.p=TRUE)
  )
  # The one thing we do need to fix is zero probability for non-integer values
  # of x
  noninteger <- recycle(x != floor(x), length(lp))
  lp[(!is.nan(lp)) & noninteger] <- -Inf
  # If we don't want the log, exponentiate lp
  if (log) {
    return(lp)
  } else {
    return(exp(lp))
  }
}

# Quantile function for the Gamma-count distribution
# Finds the smallest x such that pgc(x, lambda, alpha) >= p. Uses Vectorize
# to work with vector arguments
# The expand/search design means that it runs with O(log(x)) calls to pgc, where
# x is the final returned number
qgc <- Vectorize(
  function(p, lambda, alpha=1, lower.tail=TRUE, log.p=FALSE) {
    # Correct p for lower.tail if necessary
    if (!lower.tail) {
      if (log.p) {
        p <- logdiffexp(0, p)
      } else {
        p <- 1 - p
      }
    }
    # This variable stores the largest known value of x that is too low:
    # pgc(low.x, lambda, alpha) < p
    low.x <- -1
    # This variable stores the smallest known value of x that is large enough:
    # pgc(high.x, lambda, alpha) >= p
    high.x <- Inf
    # Expansion phase: Check values of x = 2^k - 1 until one is found that is high enough
    x <- 0
    while(pgc(x, lambda, alpha, log.p=log.p) < p) {
      low.x <- x
      x <- 2 * x + 1
    }
    high.x <- x
    # Search phase: the value is somewhere between low.x (exclusive) and high.x (inclusive)
    # high.x - low.x is guaranteed to be a power of two, so binary search is uncomplicated
    while(high.x - low.x > 1) {
      x <- (high.x + low.x) / 2
      if (pgc(x, lambda, alpha, log.p=log.p) < p) {
        low.x <- x
      } else {
        high.x <- x
      }
    }
    # low.x and high.x are one apart, meaning that high.x is the value to return
    return(high.x)
  },
  vectorize.args = c("p", "lambda", "alpha")
)

# Random number generator for gamma-count distribution
rgc <- function(n, lambda, alpha=1) {
  u <- runif(n) # Will use probability integral transform
  return(qgc(u, recycle(lambda, n), recycle(alpha, n)))
}

