# Gamma-Count with Random Start Time
# Gamma-count process where while events begin arriving at t=0, the counted
# interval of length lambda doesn't begin until a randomly chosen start time.
# This is implemented by modifying the arrival time distribution for the first
# arrival.


#' The Gamma-Count Distribution with Random Start Time
#'
#' Density, distribution function, quantile function, and random number
#' generation for the gamma-count distribution with random start time with
#' scale/mean parameter lambda and dispersion parameter alpha
#'
#' The gamma-count distribution with random start time models the count of
#' events arriving in the interval (L,L+lambda] where L is a large randomly
#' chosen positive real number and the times between successive events follow a
#' gamma(alpha, alpha) distribution. Since the expected time between events is
#' 1, lambda is approximately the mean of the distribution. This
#' distribution has the Poisson distribution as a special case when alpha = 1.
#' The random start time improves the standard gamma-count distribution by
#' removing the behavior near t=0 for under- or over-dispersed values of alpha.
#'
#' dgcrst.approx and pgcrst.approx use an approximation to the integral which
#' is relatively fast and still fairly accurate. The largest time savings is
#' with vector input.
#'
#' @param x vector of (non-negative integer) quantiles
#' @param q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param lambda vector of (positive) scale/mean parameters
#' @param alpha vector of (positive) dispersion parameters
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE, probabilities are P[X <= x], otherwise P[X > x]
#' @param diagnostics logical; if TRUE, return diagnostics from integrate()
#' @param segments number of segments to use for approximate computation
#'
#' @return dgcrst and dgcrst.approx give the (log) density, pgcrst and
#'   pgcrst.approx give the (log) distribution function, qgcrst gives the
#'   quantile function, and rgcrst generates random counts.
#'
#'   Invalid lambda or alpha will result in return value NaN, with a warning.
#'
#'   The length of the result is determined by n for rgcrst, and the length of the
#'   longest numeric argument for the other functions. The numerical arguments
#'   other than n are recycled to the length of the result. Only the first
#'   elements of the logical arguments are used.
#'
#'   In dgcrst and pgcrst, if diagnostics=TRUE, the returned value is a list with
#'   components value, abs.error, subdivisions, message. All components are
#'   vectors with the same length. 'value' is the probability or density
#'   requested, and all other values come directly from any calls
#'   to integrate(). NA is put in the diagnostics for any value for which
#'   integrate() did not need to be called.
#'
#' @seealso [integrate()] for diagnostic output details, [dgc()] for the
#'   standard gamma-count distribution, and [dft()] for the first
#'   arrival time distribution after the random start time
#'
#' @name GammaCountRST
NULL


# PMF function
#' @rdname GammaCountRST
#' @export
dgcrst <- function(x, lambda, alpha=1, log=FALSE, diagnostics=FALSE) {
  # Determine length of output and recycle inputs to that length
  n <- max(length(x), length(lambda), length(alpha))
  x.fill      <- recycle(x, n)
  lambda.fill <- recycle(lambda, n)
  alpha.fill  <- recycle(alpha, n)

  # Integrand is the joint density of first arrival and observed count,
  # assuming first arrival is before lambda.
  # Transformed with u = tau/(1+tau)
  integrand <- function(u, x, lambda, alpha) {
    dft(u / (1 - u), alpha) * dgc(x - 1, lambda - u / (1 - u), alpha) / (1 - u)^2
  }

  # Calculate the probabilities
  probs <- rep(0, n)

  # If diagnostics requested, initialize those arrays
  if (diagnostics) {
    abs.error <- rep(NA, n)
    subdivisions <- rep(NA, n)
    message <- rep(NA, n)
  }

  # Start with x == 0 cases
  zero <- (x.fill == 0)
  probs[zero] <- pft(lambda.fill[zero], alpha=alpha.fill[zero], lower.tail=FALSE)

  # Now do the x > 0 cases
  for (i in 1:n) {
    if ((x.fill[i] > 0) && (x.fill[i] == floor(x.fill[i]))) {
      # If x is greater than zero, we integrate over all first arrival times
      # between 0 and lambda
      u.max <- lambda.fill[i] / (lambda.fill[i] + 1)
      integ.res <- integrate(f=integrand, lower=0, upper=u.max,
                             x=x.fill[i], lambda=lambda.fill[i], alpha=alpha.fill[i],
                             stop.on.error=FALSE)
      if (diagnostics) {
        probs[i] <- integ.res$value
        abs.error[i] <- integ.res$abs.error
        subdivisions[i] <- integ.res$subdivisions
        message[i] <- integ.res$message
      } else {
        probs[i] <- if (integ.res$message == "OK") integ.res$value else NA
      }
    }
  }

  # Now take the log if necessary and return
  if (log) {
    probs <- log(probs)
  }
  if (diagnostics) {
    return(
      list(value=probs,
           abs.error=abs.error,
           subdivisions=subdivisions,
           message=message)
    )
  } else {
    return(probs)
  }
}

# CDF function
#' @rdname GammaCountRST
#' @export
pgcrst <- function(x, lambda, alpha=1, lower.tail=TRUE, log.p=FALSE, diagnostics=FALSE) {
  # Determine length of output and recycle inputs to that length
  n <- max(length(x), length(lambda), length(alpha))
  x.fill      <- recycle(floor(x), n)
  lambda.fill <- recycle(lambda, n)
  alpha.fill  <- recycle(alpha, n)
  # Integrand: integrating this over 0 < tau < lambda gives ccdf of x >= 1
  # i.e. P(X > x) for x >= 1. It has been transformed by u = tau / (tau + 1)
  integrand <- function(u, x, lambda, alpha) {
    dft(u / (1 - u), alpha) * pgc(x - 1, lambda - u / (1 - u), alpha, lower.tail=FALSE) / (1 - u)^2
  }

  # Calculate the probabilities. Will fill with complementary cdf, then invert
  # later if necessary
  probs <- rep(1, n)

  # If diagnostics requested, initialize the arrays to store output from integrate()
  if (diagnostics) {
    abs.error <- rep(NA, n)
    subdivisions <- rep(NA, n)
    message <- rep(NA, n)
  }

  # Start with x == 0 cases. P(X > 0) = P(tau <= lambda)
  zero <- (x.fill == 0)
  probs[zero] <- pft(lambda.fill[zero], alpha=alpha.fill[zero], lower.tail=TRUE)

  # Now fill in the nonzero cases
  for (i in 1:n) {
    if (x.fill[i] > 0) {
      # If x > 0, find ccdf by integrating integrand
      u.max <- lambda.fill[i] / (lambda.fill[i] + 1)
      integ.res <- integrate(f=integrand, lower=0, upper=u.max,
                             x=x.fill[i], lambda=lambda.fill[i], alpha=alpha.fill[i],
                             stop.on.error=FALSE)
      if (diagnostics) {
        probs[i] <- integ.res$value
        abs.error[i] <- integ.res$abs.error
        subdivisions[i] <- integ.res$subdivisions
        message[i] <- integ.res$message
      } else {
        probs[i] <- if (integ.res$message == "OK") integ.res$value else NA
      }
    }
  }

  # Correct to lower tail and/or log scale, if needed
  if (lower.tail && log.p) {
    probs <- log1p(-probs)
  } else if (lower.tail) {
    probs <- 1 - probs
  } else if (log.p) {
    probs <- log(probs)
  }
  if (diagnostics) {
    return(
      list(value=probs,
           abs.error=abs.error,
           subdivisions=subdivisions,
           message=message)
    )
  } else {
    return(probs)
  }
}

# Quantile function for the Gamma-count random start time distribution
# Finds the smallest x such that pgc(x, lambda, alpha) >= p. Uses Vectorize
# to work with vector arguments
# The expand/search design means that it runs with O(log(x)) calls to pgcrst,
# where x is the final returned number
#' @rdname GammaCountRST
#' @export
qgcrst <- Vectorize(
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
    # pgcrst(low.x, lambda, alpha) < p
    low.x <- -1
    # This variable stores the smallest known value of x that is large enough:
    # pgcrst(high.x, lambda, alpha) >= p
    high.x <- Inf
    # Expansion phase: Check values of x = 2^k - 1 until one is found that is high enough
    x <- 0
    while(pgcrst(x, lambda, alpha, log.p=log.p) < p) {
      low.x <- x
      x <- 2 * x + 1
    }
    high.x <- x
    # Search phase: the value is somewhere between low.x (exclusive) and high.x (inclusive)
    # high.x - low.x is guaranteed to be a power of two, so binary search is uncomplicated
    while(high.x - low.x > 1) {
      x <- (high.x + low.x) / 2
      if (pgcrst(x, lambda, alpha, log.p=log.p) < p) {
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

# Random number generator for gamma-count random start time distribution
#' @rdname GammaCountRST
#' @export
rgcrst <- function(n, lambda, alpha=1) {
  # Pass this off to the joint distribution rng, avoid repeated code
  rgcrst.joint(n, lambda, alpha)[,"x"]
}

