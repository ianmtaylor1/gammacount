

#' @export
darrival <- function(x, num=1, alpha=1, log=FALSE) {

  len <- max(length(x), length(num), length(alpha))
  x <- recycle(x, len)
  num <- recycle(num, len)
  alpha <- recycle(alpha, len)

  # Enumerate the special cases to fix after the fact
  resnan <- (alpha <= 0) | (num != floor(num)) | (num < 1)
  nosupport <- (!resnan) & (x < 0)

  # Calculate the mode to use as a switchover between representations
  fmode <- exp((lgamma(num * alpha) - lgamma((num - 1) * alpha)) / alpha) / alpha

  # Calculation categories
  first <- (num == 1)             # First arrival - simpler form of density
  usep <- (!first) & (x < fmode)  # Left of the mode, use lower.tail=TRUE in pgamma
  useq <- (!first) & (!usep)      # Right of the mode, use lower.tail=FALSE in pgamma

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

  # This calculates the log complementary cdf and transforms it as necessary
  # in the end to lower tail and/or nonlog version

  len <- max(length(x), length(num), length(alpha))
  x <- recycle(x, len)
  num <- recycle(num, len)
  alpha <- recycle(alpha, len)

  logccdf <- logdiffexp(
    log(num) + pgamma(x, num * alpha + 1, alpha, lower.tail=FALSE, log.p=TRUE),
    log(x)   + pgamma(x, num * alpha,     alpha, lower.tail=FALSE, log.p=TRUE)
  )

  first <- (num == 1)
  logccdf[!first] <- logdiffexp(
    logccdf[!first],
    logdiffexp(
      log(num[!first] - 1) + pgamma(x[!first], (num[!first] - 1) * alpha[!first] + 1, alpha[!first], lower.tail=FALSE, log.p=TRUE),
      log(x[!first])       + pgamma(x[!first], (num[!first] - 1) * alpha[!first],     alpha[!first], lower.tail=FALSE, log.p=TRUE)
    )
  )

  # Correct for numerical errors early when ccdf is highest
  logccdf <- pmin(logccdf, 0)

  # Correct to lower tail and exponentiate if necessary
  if (lower.tail && !log.p) {
    exp(logdiffexp(0, logccdf))
  } else if (lower.tail) {
    logdiffexp(0, logccdf)
  } else if (!log.p) {
    exp(logccdf)
  } else {
    logccdf
  }
}

#' @export
rarrival <- function(n, num=1, alpha=1) {
  num <- recycle(num, n)
  alpha <- recycle(alpha, n)

  # First time
  tau1 <- rft(n, alpha)

  # Remaining inter-arrival times
  delta <- 0
  first <- (num == 1)
  delta[!first] <- rgamma(sum(!first), (num[!first] - 1) * alpha[!first], alpha[!first])

  # Resulting arrival times
  tau1 + delta
}
