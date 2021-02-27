# Exact computation of dgcrst and pgcrst from recently-derived cdf formula

#' @export
pgcrst.exact <- function(x, lambda, alpha=1, lower.tail=TRUE, log.p=FALSE) {
  # Make all variables the right length
  n <- max(length(x), length(lambda), length(alpha))
  x <- recycle(floor(x), n)
  lambda <- recycle(lambda, n)
  alpha <- recycle(alpha, n)

  # Will store the result
  log.probs <- rep(-Inf, n)

  zero <- (x == 0)
  pos <- (x > 0)

  # P(X <= 0) = P(X = 0) = P(tau > lambda)
  log.probs[zero] <- pft(lambda[zero], alpha[zero], lower.tail=FALSE, log.p=TRUE)

  # For positive x's, use new formula
  # log.probs[pos] <- logdiffexp(
  #   logsumexp(
  #     log1p(x[pos]) + pgamma(lambda[pos], (x[pos] + 1) * alpha[pos] + 1, alpha[pos], lower.tail=FALSE, log.p=TRUE),
  #     log(lambda[pos]) + pgamma(lambda[pos], x[pos] * alpha[pos], alpha[pos], lower.tail=FALSE, log.p=TRUE)
  #   ),
  #   logsumexp(
  #     log(x[pos]) + pgamma(lambda[pos], x[pos] * alpha[pos] + 1, alpha[pos], lower.tail=FALSE, log.p=TRUE),
  #     log(lambda[pos]) + pgamma(lambda[pos], (x[pos] + 1) * alpha[pos], alpha[pos], lower.tail=FALSE, log.p=TRUE)
  #   )
  # )
  log.probs[pos] <- logdiffexp(
    logdiffexp(
      log1p(x[pos]) + pgamma(lambda[pos], (x[pos] + 1) * alpha[pos] + 1, alpha[pos], lower.tail=FALSE, log.p=TRUE),
      log(x[pos]) + pgamma(lambda[pos], x[pos] * alpha[pos] + 1, alpha[pos], lower.tail=FALSE, log.p=TRUE)
    ),
    logdiffexp(
      pgamma(lambda[pos], (x[pos] + 1) * alpha[pos], alpha[pos], lower.tail=FALSE, log.p=TRUE),
      pgamma(lambda[pos], x[pos] * alpha[pos], alpha[pos], lower.tail=FALSE, log.p=TRUE)
    ) + log(lambda[pos])
  )

  # Transform probabilities to desired form and return
  if (!lower.tail) {
    log.probs <- logdiffexp(rep(0, length(log.probs)), log.probs)
  }
  if (log.p) {
    return(log.probs)
  } else {
    return(exp(log.probs))
  }
}

#' @export
dgcrst.exact <- function(x, lambda, alpha=1, log=FALSE) {
  log.probs <- logdiffexp(
    pgcrst.exact(x, lambda, alpha, log.p=TRUE),
    pgcrst.exact(x-1, lambda, alpha, log.p=TRUE)
  )
  if (log) {
    return(log.probs)
  } else {
    return(exp(log.probs))
  }
}
