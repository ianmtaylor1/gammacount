# This file contains functions for the Gamma count distribution
# Like poisson, but arrivals are gamma(alpha,beta). (expected value = alpha/beta)
# Length of count interval is lambda, like poisson.
# poisson is a special case with alpha=1.

#' The Gamma-Count Distribution
#'
#' Density, distribution function, quantile function, and random number
#' generation for the gamma-count distribution with time parameter lambda
#' and gamma distribution parameters alpha and beta.
#'
#' The gamma-count distribution models the count of events arriving in the
#' interval (0,lambda] when the times between successive events follow a
#' gamma(alpha, beta) distribution. This is a renewal process with gamma renewal
#' times. lambda * beta / alpha is approximately the mean. This distribution
#' has the Poisson distribution as a special case when alpha = beta = 1.
#'
#' Note that lambda and beta have the same effect on the count distribution. The
#' distribution of counts only depends on either through the product
#' lambda * beta, so only one needs to be set. Therefore by default,
#' beta = alpha. This way arrival times have mean 1 and lambda can be used to
#' approximately model the mean of the count distribution.
#'
#' @param x vector of (non-negative integer) quantiles
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param lambda vector of (positive) time length parameters
#' @param alpha vector of (positive) dispersion parameters
#' @param beta vector of (positive) rate parameters. By default, beta = alpha
#'   so that arrival times have mean 1.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE, probabilities are P[X <= x], otherwise P[X > x]
#'
#' @return dgc gives the (log) density, pgc gives the (log) distribution
#'   function, qgc gives the quantile function, and rgc generates random counts.
#'
#'   Invalid lambda, alpha, or beta will result in return value NaN.
#'
#'   The length of the result is determined by n for rgc, and the length of the
#'   longest numeric argument for the other functions. The numerical arguments
#'   other than n are recycled to the length of the result. Only the first
#'   elements of the logical arguments are used.
#'
#' @name GammaCount
NULL

# Probability mass funtion for gamma-count distribution
#' @rdname GammaCount
#' @export
dgc <- function(x, lambda, alpha=1, beta=alpha, log=FALSE) {
  n <- max(length(x), length(lambda), length(alpha), length(beta))
  x <- recycle(x, n)
  lambda <- recycle(lambda, n)
  alpha <- recycle(alpha, n)
  beta <- recycle(beta, n)
  # Define the pmf in terms of the cdf. Choose direction based on relation of x
  # to the approximate center of the distribution.
  left <- (x <= lambda * beta / alpha)
  left[is.na(left)] <- FALSE
  lp <- rep(NaN, n)
  lp[left] <- logdiffexp(
    pgc(x[left],     lambda[left], alpha[left], beta[left], log.p=TRUE),
    pgc(x[left] - 1, lambda[left], alpha[left], beta[left], log.p=TRUE)
  )
  lp[!left] <- logdiffexp(
    pgc(x[!left] - 1, lambda[!left], alpha[!left], beta[!left], log.p=TRUE, lower.tail=FALSE),
    pgc(x[!left],     lambda[!left], alpha[!left], beta[!left], log.p=TRUE, lower.tail=FALSE)
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

# Cumulative distribution function for gamma-count distribution
#' @rdname GammaCount
#' @export
pgc <- function(q, lambda, alpha=1, beta=alpha, lower.tail=TRUE, log.p=FALSE) {
  # Use pgamma to calculate the probabilities (following Zeviani et al. (2013))
  # A few useful tricks:
  #  1. Fill invalid parameter values with NaN so that the NaN propogates to the
  #     result.
  #  2. pgamma with a zero shape parameter is kind of a point mass at zero, so
  #     x == -1 behaves nicely for all x < 0. Use pmax to provide lower bound.
  pgamma(ifelse(lambda > 0, lambda, NaN),
         (pmax(floor(q), -1) + 1) * ifelse(alpha > 0, alpha, NaN),
         rate = ifelse(beta > 0, beta, NaN),
         lower.tail = !lower.tail,
         log.p = log.p)
}

# Quantile function for the Gamma-count distribution
# Finds the smallest x such that pgc(x, lambda, alpha) >= p. Uses Vectorize
# to work with vector arguments
# The expand/search design means that it runs with O(log(x)) calls to pgc, where
# x is the final returned number
#' @rdname GammaCount
#' @export
qgc <- Vectorize(
  function(p, lambda, alpha=1, beta=alpha, lower.tail=TRUE, log.p=FALSE) {
    # Check for high and low values of p to return Inf or NaN
    if (lower.tail && !log.p) {
      if (p == 1) return(Inf)
      if ((p > 1) || (p < 0)) return(NaN)
    }
    if (lower.tail && log.p) {
      if (p == 0) return(Inf)
      if (p > 0) return(NaN)
    }
    if (!lower.tail && !log.p) {
      if (p == 0) return(Inf)
      if ((p > 1) || (p < 0)) return(NaN)
    }
    if (!lower.tail && log.p) {
      if (p == -Inf) return(Inf)
      if (p > 0) return (NaN)
    }
    # Check for invalid values of parameters to return NaN
    if ((lambda <= 0) || (alpha <= 0) || (beta <= 0)) return(NaN)

    # Correct p for lower.tail if necessary. Value multiplied by p and pgc().
    mult <- if (lower.tail) { 1 } else { -1 }

    # low.x stores the largest known value of x that is too low:
    # pgc(low.x, lambda, alpha) < p
    # high.x stores the smallest known value of x that is large enough:
    # pgc(high.x, lambda, alpha) >= p

    # Expansion phase: search out from zero until finding a
    # range such that q \in (low.x, high.x] and the two above inequalities are
    # satisfied. Only one of expand-up or expand-down will ever happen
    low.x <- -1
    high.x <- 0
    while(pgc(high.x, lambda, alpha, beta, log.p=log.p, lower.tail=lower.tail) * mult < p * mult) {
      low.x <- high.x
      high.x <- high.x * 2 + 1
    }

    # Search phase: the value is somewhere between low.x (exclusive) and high.x
    # (inclusive). Divide this range roughly in half until the value is found
    while(high.x - low.x > 1) {
      x <- (high.x + low.x) %/% 2
      if (pgc(x, lambda, alpha, beta, log.p=log.p, lower.tail=lower.tail) * mult < p * mult) {
        low.x <- x
      } else {
        high.x <- x
      }
    }

    # low.x and high.x are one apart, meaning that high.x is the value to return
    return(high.x)
  },
  vectorize.args = c("p", "lambda", "alpha", "beta")
)

# Random number generator for gamma-count distribution
#' @rdname GammaCount
#' @export
rgc <- function(n, lambda, alpha=1, beta=alpha) {
  u <- runif(n) # Will use probability integral transform
  return(qgc(u, recycle(lambda, n), recycle(alpha, n), recycle(beta, n)))
}

