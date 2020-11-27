# This file contains functions to approximate the more computationally expensive
# Gammacount Random Start Time distribuiton functions. All end in 'a' denoting
# their approximate nature


# Determine the bounds of segments used to approximate the integral
# Take a value p between 0 and 1, and return the corresponding value between
# 0 and lambda for creating segments.
# Works somewhat like the following:
#  Divide the interval [0,1] evenly between segments. Divide the remaining
#  length, lambda-1 assuming lambda > 1, logarithmically with base = exp(alpha)
find.bound <- function(p, lambda, alpha) {
  max.p <- -expm1(-alpha * pmax(lambda - 1, 0))
  return(
    pmin(p * lambda, p - log1p(-p * max.p) / alpha)
  )
}


# Approximate dgcrst function that doesn't have to call integrate()
#' @rdname GammaCountRST
#' @export
dgcrst.approx <- function(x, lambda, alpha=1, log=FALSE, segments=10) {
  # Determine length of output and recycle inputs to that length
  n <- max(length(x), length(lambda), length(alpha))
  x.fill      <- recycle(x, n)
  lambda.fill <- recycle(lambda, n)
  alpha.fill  <- recycle(alpha, n)

  # Calculate the probabilities
  log.probs <- rep(-Inf, n)

  # Start with x == 0 cases
  zero <- (x.fill == 0)
  pos <- (x.fill > 0) & (x.fill == floor(x.fill))
  log.probs[zero] <- pft(lambda.fill[zero], alpha.fill[zero], lower.tail=FALSE, log.p=TRUE)

  # Now do all the positive x's
  # Divide the interval [0, lambda] into 'segments' number of segments spaced
  # depending on alpha. Hihger alpha (more underdispersed) and more of the
  # segments will be concentrated between 0 and 1, lower alpha (more overdispersed)
  # and the segments will be more evenly spread between 0 and lambda.
  # Then for each segment, assume the first arrival time was the conditional
  # expected value within that segment, and then the result is a weighted average
  # of dgc(..., lambda - tau, ...) based on the probability of each segment.
  log.cdf.prev <- -Inf
  log.partev.prev <- -Inf
  for (seg in seq_len(segments)/segments) {
    # What is the upper bound arrival time for this segment?
    seg.bound <- find.bound(seg, lambda.fill[pos], alpha.fill[pos])
    # What is the expected first arrival time between the last segment boundary
    # and this one?
    log.cdf <- pft(seg.bound, alpha.fill[pos], log.p=TRUE)
    log.partev <- ft.partial.ev(seg.bound, alpha.fill[pos], log=TRUE)
    tau <- exp(logdiffexp(log.partev, log.partev.prev) - logdiffexp(log.cdf, log.cdf.prev))
    # Conditioned on tau, accumulate the probability of x
    log.probs[pos] <- logsumexp(
      log.probs[pos],
      logdiffexp(log.cdf, log.cdf.prev) + dgc(x.fill[pos] - 1, lambda.fill[pos] - tau, alpha.fill[pos], log=TRUE)
    )
    # Switch the previous values of cdf and partial expected value before the
    # next loop
    log.cdf.prev <- log.cdf
    log.partev.prev <- log.partev
  }

  # Now take the log if necessary and return
  if (log) {
    return(log.probs)
  } else {
    return(exp(log.probs))
  }
}


# Approximate pgcrst function that doesn't have to call integrate()
#' @rdname GammaCountRST
#' @export
pgcrst.approx <- function(x, lambda, alpha=1, log.p=FALSE, lower.tail=TRUE, segments=10) {
  # Determine length of output and recycle inputs to that length
  n <- max(length(x), length(lambda), length(alpha))
  x.fill      <- recycle(floor(x), n)
  lambda.fill <- recycle(lambda, n)
  alpha.fill  <- recycle(alpha, n)

  # Calculate the probabilities. We will be building up the log upper tail by
  # default, then transforming later
  log.probs <- rep(0, n)

  # Start with x == 0 cases
  zero <- (x.fill == 0)
  pos <- (x.fill > 0)
  log.probs[zero] <- pft(lambda.fill[zero], alpha.fill[zero], lower.tail=TRUE, log.p=TRUE)

  # Now do all the positive x's
  # Divide the interval [0, lambda] into 'segments' number of segments spaced
  # depending on alpha. Hihger alpha (more underdispersed) and more of the
  # segments will be concentrated between 0 and 1, lower alpha (more overdispersed)
  # and the segments will be more evenly spread between 0 and lambda.
  # Then for each segment, assume the first arrival time was the conditional
  # expected value within that segment, and then the result is a weighted average
  # of pgc(..., lambda - tau, ...) based on the probability of each segment.
  log.cdf.prev <- -Inf
  log.partev.prev <- -Inf
  log.probs[pos] <- -Inf  # Log-scale sum depends on this starting value
  for (seg in seq_len(segments)/segments) {
    # What is the upper bound arrival time for this segment?
    seg.bound <- find.bound(seg, lambda.fill[pos], alpha.fill[pos])
    # What is the expected first arrival time between the last segment boundary
    # and this one?
    log.cdf <- pft(seg.bound, alpha.fill[pos], log.p=TRUE)
    log.partev <- ft.partial.ev(seg.bound, alpha.fill[pos], log=TRUE)
    tau <- exp(logdiffexp(log.partev, log.partev.prev) - logdiffexp(log.cdf, log.cdf.prev))
    # Conditioned on tau, accumulate the probability of x
    log.probs[pos] <- logsumexp(
      log.probs[pos],
      logdiffexp(log.cdf, log.cdf.prev) + pgc(x.fill[pos] - 1, lambda.fill[pos] - tau, alpha.fill[pos], log.p=TRUE, lower.tail=FALSE)
    )
    # Switch the previous values of cdf and partial expected value before the
    # next loop
    log.cdf.prev <- log.cdf
    log.partev.prev <- log.partev
  }

  # Now take the log if necessary and return
  if (lower.tail && (!log.p)) {
    return(-expm1(log.probs))
  } else if (lower.tail) {
    return(logdiffexp(0, log.probs))
  } else if (!log.p) {
    return(exp(log.probs))
  } else {
    return(log.probs)
  }
}
