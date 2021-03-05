

#' @export
darrival <- function(x, num=1, alpha=1, log=FALSE) {

  len <- max(length(x), length(num), length(alpha))
  x <- recycle(x, len)
  num <- recycle(num, len)
  alpha <- recycle(alpha, len)

  # Enumerate the special cases to fix after the fact
  resnan <- (!is.finite(alpha)) | (!is.finite(num)) | (alpha <= 0) | (num != floor(num)) | (num < 1)
  nosupport <- (!resnan) & (x < 0)

  # Calculate the mode to use as a switchover between representations
  fmode <- exp((lgamma(num * alpha) - lgamma((num - 1) * alpha)) / alpha) / alpha

  # Calculation categories
  first <- (!resnan) & (num == 1)             # First arrival - simpler form of density
  usep <- (!resnan) & (!first) & (x < fmode)  # Left of the mode, use lower.tail=TRUE in pgamma
  useq <- (!resnan) & (!first) & (!usep)      # Right of the mode, use lower.tail=FALSE in pgamma

  logdens <- rep(NaN, len)

  # For first arrivals, use only one pgamma
  logdens[first] <- pgamma(x[first], alpha[first], alpha[first], lower.tail=FALSE, log.p=TRUE)
  # For left-hand, use pgamma with lower tail
  logdens[usep] <- logdiffexp(
    pgamma(x[usep], (num[usep] - 1) * alpha[usep], alpha[usep], lower.tail=TRUE, log.p=TRUE),
    pgamma(x[usep], num[usep] * alpha[usep],       alpha[usep], lower.tail=TRUE, log.p=TRUE)
  )
  # For right-hand, use pgamma with upper tail
  logdens[useq] <- logdiffexp(
    pgamma(x[useq], num[useq] * alpha[useq],       alpha[useq], lower.tail=FALSE, log.p=TRUE),
    pgamma(x[useq], (num[useq] - 1) * alpha[useq], alpha[useq], lower.tail=FALSE, log.p=TRUE)
  )

  # Fix edge cases
  logdens[resnan] <- NaN
  logdens[nosupport] <- -Inf

  if (log) {
    logdens
  } else {
    exp(logdens)
  }
}


#' @export
parrival <- function(x, num=1, alpha=1, lower.tail=TRUE, log.p=FALSE) {

  len <- max(length(x), length(num), length(alpha))
  x <- recycle(pmax(x, 0), len)
  num <- recycle(num, len)
  alpha <- recycle(alpha, len)

  # Special case to return NaN - invalid parameter values
  resnan <- (!is.finite(alpha)) | (!is.finite(num)) | (alpha <= 0) | (num != floor(num)) | (num < 1)

  # Fill in the non-NaN values of the logcdf from the appropriate subfunction
  logcdf <- rep(NaN, len)
  first <- ((!resnan) & (num == 1))
  later <- ((!resnan) & (num > 1))
  if (lower.tail) {
    logcdf[first] <- parrival_left_first(x[first], alpha[first])
    logcdf[later] <- parrival_left_later(x[later], num[later], alpha[later])
  } else {
    logcdf[first] <- parrival_right_first(x[first], alpha[first])
    logcdf[later] <- parrival_right_later(x[later], num[later], alpha[later])
  }

  # Check if recomputing in the opposite tail might be beneficial
  # switchpoint value MUST be larger than log(1/2) = -0.693 to avoid infinite
  # recursion, and preferably leave a little room
  switchpoint <- -0.5
  recompute <- (!resnan) & (logcdf > switchpoint)
  if (sum(recompute) > 0) {
    logcdf[recompute] <- logdiffexp(
      0,
      parrival(x[recompute], num[recompute], alpha[recompute], log.p=TRUE, lower.tail=!lower.tail)
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

#' @export
rarrival <- function(n, num=1, alpha=1) {
  num <- recycle(num, n)
  alpha <- recycle(alpha, n)

  resnan <- (!is.finite(alpha)) | (!is.finite(num)) | (alpha <= 0) | (num != floor(num)) | (num < 1)

  # First time
  tau1 <- runif(n, 0, rgamma(n, alpha + 1, alpha))

  # Remaining inter-arrival times
  delta <- rep(0, n)
  first <- (num == 1)
  delta[!first] <- rgamma(sum(!first), (num[!first] - 1) * alpha[!first], alpha[!first])

  # Resulting arrival times
  val <- tau1 + delta
  val[resnan] <- NaN
  val
}


######## Sub-functions for parrival in different cases
# Each returns the log cdf as described in its top comment.
# Each does no parameter checking, they strictly do the computation.

# Left-hand probability for first arrivals (num == 1)
parrival_left_first <- function(x, alpha) {
  logsumexp(
    pgamma(x, alpha+1, rate=alpha, log.p=TRUE, lower.tail=TRUE),
    log(x) + pgamma(x, alpha, rate=alpha, log.p=TRUE, lower.tail=FALSE)
  )
}

# Left-hand probability for later arrivals (num > 1)
parrival_left_later <- function(x, num, alpha) {
  logdiffexp(
    logsumexp(
      log(num) + pgamma(x, num * alpha + 1,   alpha, log.p=TRUE, lower.tail=TRUE),
      log(x)   + pgamma(x, (num - 1) * alpha, alpha, log.p=TRUE, lower.tail=TRUE)
    ),
    logsumexp(
      log(num - 1) + pgamma(x, (num - 1) * alpha + 1, alpha, log.p=TRUE, lower.tail=TRUE),
      log(x)       + pgamma(x, num * alpha,           alpha, log.p=TRUE, lower.tail=TRUE)
    )
  )
}

# Right-hand probability for first arrivals (num == 1)
parrival_right_first <- function(x, alpha) {
  logdiffexp(
    pgamma(x, alpha+1, rate=alpha, log.p=TRUE, lower.tail=FALSE),
    log(x) + pgamma(x, alpha, rate=alpha, log.p=TRUE, lower.tail=FALSE)
  )
}

# Right-hand probability for later arrivals (num > 1)
parrival_right_later <- function(x, num, alpha) {
  logdiffexp(
    logsumexp(
      log(num) + pgamma(x, num * alpha + 1,   alpha, log.p=TRUE, lower.tail=FALSE),
      log(x)   + pgamma(x, (num - 1) * alpha, alpha, log.p=TRUE, lower.tail=FALSE)
    ),
    logsumexp(
      log(num - 1) + pgamma(x, (num - 1) * alpha + 1, alpha, log.p=TRUE, lower.tail=FALSE),
      log(x)       + pgamma(x, num * alpha,           alpha, log.p=TRUE, lower.tail=FALSE)
    )
  )
}
