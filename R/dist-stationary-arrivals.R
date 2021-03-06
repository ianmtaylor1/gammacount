#' Stationary Gamma Arrival distribution
#'
#' Density, distribution function, quantile function, and random number
#' generation for the arrival/renewal times in a stationary gamma renewal
#' process with parameters alpha and beta.
#'
#' The stationary gamma renewal process is defined as follows: the first renewal
#' time is distributed according to the limiting excess life distribution in a
#' gamma renewal process. The times between all subsequent renewal times are
#' independently distributed gamma(alpha, beta).
#'
#' These functions define the distribution of these arrival or renewal times,
#' the first of which is the limiting excess life distribution.
#'
#' The 'glel' functions are shorthand for the 'sga' functions with num = 1.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param num vector of positive integers for the number of the arrival in
#'   sequence. 1 is first, 2 is second, and so on.
#' @param alpha vector of (positive) dispersion parameters
#' @param beta vector of (positive) rate parameters. By default, beta = alpha
#'   so that arrival times have mean 1.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE, probabilities are P[X <= x], otherwise P[X > x]
#' @param tol the tolerance for calculating quantiles
#' @param max.steps the maximum number of steps used to compute the quantiles
#'
#' @return dsga and dglel give the (log) density, psga and plel give the
#'   (log) distribution function, qsga and qlel give the quantile function, and
#'   rsga and rlel generate random times.
#'
#'   Invalid alpha or beta will result in return value NaN.
#'
#'   The length of the result is determined by n for rsga and rlel, and the
#'   length of the longest numeric argument for the other functions. The
#'   numerical arguments other than n are recycled to the length of the result.
#'   Only the first elements of the logical arguments are used.
#'
#'   The returned values by qsga and qlel are approximate, but are guaranteed
#'   to be within tol of the true value. A bisection search is used, so
#'   computation time increases with -log(1/tol). If max.steps is reached, the
#'   current estimate is returned with a warning.
#'
#' @name StationaryGammaArrival
NULL

#' @rdname StationaryGammaArrival
#' @export
dsga <- function(x, num=1, alpha=1, beta=alpha, log=FALSE) {

  len <- max(length(x), length(num), length(alpha), length(beta))
  x <- recycle(x, len)
  num <- recycle(num, len)
  alpha <- recycle(alpha, len)
  beta <- recycle(beta, len)

  # Enumerate the special cases to fix after the fact
  resnan <- (
    (!is.finite(alpha))  | (alpha <= 0)
    | (!is.finite(num))  | (num != floor(num)) | (num < 1)
    | (!is.finite(beta)) | (beta <= 0)
  )
  nosupport <- (!resnan) & (x < 0)

  # Calculate the mode to use as a switchover between representations
  fmode <- exp((lgamma(num * alpha) - lgamma((num - 1) * alpha)) / alpha) / beta

  # Calculation categories
  first <- (!resnan) & (num == 1)             # First arrival - simpler form of density
  usep <- (!resnan) & (!first) & (x < fmode)  # Left of the mode, use lower.tail=TRUE in pgamma
  useq <- (!resnan) & (!first) & (!usep)      # Right of the mode, use lower.tail=FALSE in pgamma

  logdens <- rep(NaN, len)

  # For first arrivals, use only one pgamma
  logdens[first] <- pgamma(x[first], alpha[first], beta[first], lower.tail=FALSE, log.p=TRUE) - log(alpha / beta)
  # For left-hand, use pgamma with lower tail
  logdens[usep] <- logdiffexp(
    pgamma(x[usep], (num[usep] - 1) * alpha[usep], beta[usep], lower.tail=TRUE, log.p=TRUE),
    pgamma(x[usep], num[usep] * alpha[usep],       beta[usep], lower.tail=TRUE, log.p=TRUE)
  ) - log(alpha / beta)
  # For right-hand, use pgamma with upper tail
  logdens[useq] <- logdiffexp(
    pgamma(x[useq], num[useq] * alpha[useq],       alpha[useq], lower.tail=FALSE, log.p=TRUE),
    pgamma(x[useq], (num[useq] - 1) * alpha[useq], alpha[useq], lower.tail=FALSE, log.p=TRUE)
  ) - log(alpha / beta)

  # Fix edge cases
  logdens[resnan] <- NaN
  logdens[nosupport] <- -Inf

  if (log) {
    logdens
  } else {
    exp(logdens)
  }
}

#' @rdname StationaryGammaArrival
#' @export
psga <- function(x, num=1, alpha=1, beta=alpha, lower.tail=TRUE, log.p=FALSE) {

  len <- max(length(x), length(num), length(alpha), length(beta))
  x <- recycle(pmax(x, 0), len)
  num <- recycle(num, len)
  alpha <- recycle(alpha, len)
  beta <- recycle(beta, len)

  # Special case to return NaN - invalid parameter values
  resnan <- (
    (!is.finite(alpha))  | (alpha <= 0)
    | (!is.finite(num))  | (num != floor(num)) | (num < 1)
    | (!is.finite(beta)) | (beta <= 0)
  )

  # Fill in the non-NaN values of the logcdf from the appropriate subfunction
  logcdf <- rep(NaN, len)
  first <- ((!resnan) & (num == 1))
  later <- ((!resnan) & (num > 1))
  if (lower.tail) {
    logcdf[first] <- psga_left_first(x[first], alpha[first], beta[first])
    logcdf[later] <- psga_left_later(x[later], num[later], alpha[later], beta[later])
  } else {
    logcdf[first] <- psga_right_first(x[first], alpha[first], beta[first])
    logcdf[later] <- psga_right_later(x[later], num[later], alpha[later], beta[later])
  }

  # Check if recomputing in the opposite tail might be beneficial
  # switchpoint value MUST be larger than log(1/2) = -0.693 to avoid infinite
  # recursion, and preferably leave a little room
  switchpoint <- -0.5
  recompute <- (!resnan) & (logcdf > switchpoint)
  if (sum(recompute) > 0) {
    logcdf[recompute] <- logdiffexp(
      0,
      psga(x[recompute], num[recompute], alpha[recompute], beta[recompute], log.p=TRUE, lower.tail=!lower.tail)
    )
  }

  # Correct for numerical errors, if present
  logcdf <- pmin(logcdf, 0)

  # Exponentiate if necessary and return
  if (log.p) {
    logcdf
  } else {
    exp(logcdf)
  }
}

#' @rdname StationaryGammaArrival
#' @export
qsga <- function(p, num=1, alpha=1, beta=alpha, lower.tail=TRUE, log.p=FALSE,
                 tol=.Machine$double.eps ^ 0.25, max.steps=1000) {
  # Fill the arguments like other q___ functions do
  len <- max(length(p), length(num), length(alpha), length(beta))
  p <- recycle(p, len)
  num <- recycle(num, len)
  alpha <- recycle(alpha, len)
  beta <- recycle(beta, len)

  # Special case to return NaN - invalid parameter values
  resnan <- (
    (!is.finite(alpha))  | (alpha <= 0)
    | (!is.finite(num))  | (num != floor(num)) | (num < 1)
    | (!is.finite(beta)) | (beta <= 0)
  )

  # Setup: bounds and convergence tracking
  # This distribution is stochastically greater than gamma((num-1)*alpha, beta)
  # because it equals that plus a positive r.v.
  x.min <- qgamma(p, (num - 1) * alpha, rate=beta, log.p=log.p, lower.tail=lower.tail)
  # This distribution is stochactically less than gamma(num*a+1, beta) because
  # the limiting distribution of excess life is stochastically less than
  # gamma(alpha+1, beta).
  x.max <- qgamma(p, num * alpha + 1, rate=beta, log.p=log.p, lower.tail=lower.tail)

  # The previous call to qgamma checks p and alpha for us. If p is 1 or outside
  # the range (0,1) x.max will be Inf or NaN, respectively. If alpha < 0, then
  # x.max will be NaN. Either value will carry through in the final step of
  # (x.min + x.max) / 2
  done <- resnan | (!is.finite(x.max)) | (x.max - x.min <= 2 * tol)

  # Correct p for lower.tail if necessary. Value multiplied by p and pft().
  mult <- if (lower.tail) { 1 } else { -1 }

  # Loop until everyone has converged
  step <- 1
  while ((!all(done)) && (step <= max.steps)) {
    x <- (x.min + x.max) / 2
    # fx is a vector of evaluated function values with length == sum(!done)
    # We avoid recalling pft again for entries which are already sufficiently
    # converged
    fx <- (psga(x[!done], num[!done], alpha[!done], beta[!done], log.p=log.p, lower.tail=lower.tail) - p[!done]) * mult
    # If fx <= 0, then x is too low of a point
    x.min[!done][fx <= 0] <- x[!done][fx <= 0]
    # If fx >= 0, then x is too high of a point
    x.max[!done][fx >= 0] <- x[!done][fx >= 0]
    # See who is done now
    done <- done | (x.max - x.min <= 2 * tol) # Things only get more done.
    step <- step + 1
  }
  if ((step > max.steps) && (!all(done))) {
    warning("Maximum steps exceeded, some results may not be within tolerance.")
  }
  # By returning the average of x.min and x.max (which are at most 2 * tol
  # apart), we guarantee the result is within tol of a root
  return((x.min + x.max) / 2)
}

#' @rdname StationaryGammaArrival
#' @export
rsga <- function(n, num=1, alpha=1, beta=alpha) {
  num <- recycle(num, n)
  alpha <- recycle(alpha, n)
  beta <- recycle(beta, n)

  resnan <- (
    (!is.finite(alpha))  | (alpha <= 0)
    | (!is.finite(num))  | (num != floor(num)) | (num < 1)
    | (!is.finite(beta)) | (beta <= 0)
  )

  # First time
  tau1 <- runif(n, 0, rgamma(n, alpha + 1, beta))

  # Remaining inter-arrival times
  delta <- rep(0, n)
  first <- (num == 1)
  delta[!first] <- rgamma(sum(!first), (num[!first] - 1) * alpha[!first], beta[!first])

  # Resulting arrival times
  val <- tau1 + delta
  val[resnan] <- NaN
  val
}


################################################################################
# Gamma Limiting Excess Life distribution
#
# The distribution of excess life for a gamma renewal process as time goes to
# infinity. This is also the distribution of the first arrival/renewal time in
# a stationary gamma renewal process.

#' @rdname StationaryGammaArrival
#' @export
dglel <- function(x, alpha=1, beta=alpha, log=FALSE) {
  dsga(x, 1, alpha, beta, log)
}

#' @rdname StationaryGammaArrival
#' @export
pglel <- function(x, alpha=1, beta=alpha, lower.tail=TRUE, log.p=FALSE) {
  psga(x, 1, alpha, beta, lower.tail=lower.tail, log.p=log.p)
}

#' @rdname StationaryGammaArrival
#' @export
qglel <- function(p, alpha=1, beta=alpha, lower.tail=TRUE, log.p=FALSE,
                  tol=.Machine$double.eps ^ 0.25, max.steps=1000) {
  qsga(p, 1, alpha, beta, lower.tail=lower.tail, log.p=log.p,
       tol=tol, max.steps=max.steps)
}

#' @rdname StationaryGammaArrival
#' @export
rglel <- function(n, alpha=1, beta=alpha) {
  rsga(n, 1, alpha, beta)
}

################################################################################

######## Sub-functions for psga in different cases
# Each returns the log cdf as described in its top comment.
# Each does no parameter checking, they strictly do the computation.

# Left-hand probability for first arrivals (num == 1)
psga_left_first <- function(x, alpha, beta) {
  logsumexp(
    pgamma(x, alpha + 1, rate=beta, log.p=TRUE, lower.tail=TRUE),
    log(x * (beta / alpha)) + pgamma(x, alpha, rate=beta, log.p=TRUE, lower.tail=FALSE)
  )
}

# Left-hand probability for later arrivals (num > 1)
psga_left_later <- function(x, num, alpha, beta) {
  logdiffexp(
    logsumexp(
      log(num)                + pgamma(x, num * alpha + 1,   beta, log.p=TRUE, lower.tail=TRUE),
      log(x * (beta / alpha)) + pgamma(x, (num - 1) * alpha, beta, log.p=TRUE, lower.tail=TRUE)
    ),
    logsumexp(
      log(num - 1)            + pgamma(x, (num - 1) * alpha + 1, beta, log.p=TRUE, lower.tail=TRUE),
      log(x * (beta / alpha)) + pgamma(x, num * alpha,           beta, log.p=TRUE, lower.tail=TRUE)
    )
  )
}

# Right-hand probability for first arrivals (num == 1)
psga_right_first <- function(x, alpha, beta) {
  logdiffexp(
    pgamma(x, alpha + 1, rate=beta, log.p=TRUE, lower.tail=FALSE),
    log(x * (beta / alpha)) + pgamma(x, alpha, rate=beta, log.p=TRUE, lower.tail=FALSE)
  )
}

# Right-hand probability for later arrivals (num > 1)
psga_right_later <- function(x, num, alpha, beta) {
  logdiffexp(
    logsumexp(
      log(num)                + pgamma(x, num * alpha + 1,   beta, log.p=TRUE, lower.tail=FALSE),
      log(x * (beta / alpha)) + pgamma(x, (num - 1) * alpha, beta, log.p=TRUE, lower.tail=FALSE)
    ),
    logsumexp(
      log(num - 1)            + pgamma(x, (num - 1) * alpha + 1, beta, log.p=TRUE, lower.tail=FALSE),
      log(x * (beta / alpha)) + pgamma(x, num * alpha,           beta, log.p=TRUE, lower.tail=FALSE)
    )
  )
}
