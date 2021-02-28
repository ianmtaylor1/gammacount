

#' @export
darrival <- function(x, num=1, alpha=1, log=FALSE) {

  len <- max(length(x), length(num), length(alpha))
  x <- recycle(x, len)
  num <- recycle(num, len)
  alpha <- recycle(alpha, len)

  logdens <- pgamma(x, num * alpha, alpha, lower.tail=FALSE, log.p=TRUE)

  first <- (num == 1)
  logdens[!first] <- logdiffexp(
    logdens[!first],
    pgamma(x[!first], (num[!first] - 1) * alpha[!first], alpha[!first], lower.tail=FALSE, log.p=TRUE)
  )

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
