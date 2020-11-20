# First arrival Time distribution
# Distribution of first arrival time of a gamma-count process with random
# observation start time L >> 0

# pdf of the first time distribution
dft <- function(x, alpha=1, log=FALSE) {
  pgamma(x, shape=alpha, rate=alpha, lower.tail=FALSE, log.p=log)
}

# CDF function
pft <- function(x, alpha=1, lower.tail=TRUE, log.p=FALSE) {
  # Cap x on the low end at zero.
  x <- pmax(x, 0)
  # Calculate the log of either the lower or upper tail
  # See derivation in statdocs
  lp <- if (lower.tail) {
    logsumexp(
      pgamma(x, alpha+1, rate=alpha, log.p=TRUE, lower.tail=TRUE),
      log(x) + pgamma(x, alpha, rate=alpha, log.p=TRUE, lower.tail=FALSE)
    )
  } else {
    logdiffexp(
      pgamma(x, alpha+1, rate=alpha, log.p=TRUE, lower.tail=FALSE),
      log(x) + pgamma(x, alpha, rate=alpha, log.p=TRUE, lower.tail=FALSE)
    )
  }
  # Return the probability on the appropriate scale
  if (log.p) {
    return(lp)
  } else {
    return(exp(lp))
  }
}

# Random number generator for first arrival time distribution in gamma-count
# process
rft <- function(n, alpha=1) {
  # First draw the length of the gap we landed in
  # alpha + 1 shape because longer gaps are more likely, proportional to length
  gaps <- rgamma(n, shape=alpha + 1, rate=alpha)
  # Assume uniform within that gap
  runif(n, min=0, max=gaps)
}

# Quantile function for first arrival time distribution, using a
# faster method by using a home-spun binary search than the old optimze()
# way of doing it.
qft <- function(p, alpha=1, lower.tail=TRUE, log.p=FALSE,
                 tol=.Machine$double.eps ^ 0.25, max.steps=1000) {
  # Correct p to put it on the lower tail, if necessary
  if (!lower.tail) {
    if (log.p) {
      p <- logdiffexp(0, p)
    } else {
      p <- 1 - p
    }
  }
  # Fill the arguments like other q___ functions do
  n <- max(length(p), length(alpha))
  p <- recycle(p, n)
  alpha <- recycle(alpha, n)

  # Setup: bounds and convergence tracking
  # This distribution is stochastically greater than unif(0,1) (its density
  # starts at one and is decreasing, so an integral from 0 to x of both
  # the ft density and the uniform density will always be greater for the
  # uniform density).
  x.min <- if (log.p) { exp(p) } else { p }
  # This distribution is stochactically less than gamma(a+1, a) due to how
  # it is defined. So the p quantile of gamma(a+1, a) will always be greater
  # than the p quantile of this distribution
  x.max <- qgamma(p, alpha+1, rate=alpha, log.p=log.p)
  done <- (x.max - x.min <= 2 * tol)

  # Loop until everyone has converged
  step <- 1
  while ((!all(done)) && (step <= max.steps)) {
    x <- (x.min + x.max) / 2
    # fx is a vector of evaluated function values with length == sum(!done)
    # We avoid recalling pft again for entries which are already sufficiently
    # converged
    fx <- pft(x[!done], alpha[!done], log.p=log.p) - p[!done]
    # If fx <= 0, then x is too low of a point
    x.min[!done][fx <= 0] <- x[!done][fx <= 0]
    # If fx >= 0, then x is too high of a point
    x.max[!done][fx >= 0] <- x[!done][fx >= 0]
    # See who is done now
    done <- (x.max - x.min <= 2 * tol)
    step <- step + 1
  }
  if ((step > max.steps) && (!all(done))) {
    warning("Maximum steps exceeded, some results may not be within tolerance.")
  }
  # By returning the average of x.min and x.max (which are at most 2 * tol
  # apart), we guarantee the result is within tol of a root
  return((x.min + x.max) / 2)
}

# Calculate the "partial expected value" of the first arrival time distribution.
# Calculates P(tau <= x) * E[tau | tau <= x] (for lower tail, or greater than x
# for upper tail). This is equivalent to \int_0^x \tau f(\tau) d\tau for
# lower tail, or with bounds from x to \infty for upper tail.
ft.partial.ev <- function(x, alpha=1, lower.tail=TRUE, log=FALSE) {
  # See derivation in statdocs
  log.pev <- if (lower.tail) {
    logsumexp(
      log1p(alpha) - log(alpha) + pgamma(x, alpha + 2, rate=alpha, log.p=TRUE, lower.tail=TRUE),
      2 * log(x) + pgamma(x, alpha, rate=alpha, log.p=TRUE, lower.tail=FALSE)
    ) - log(2)
  } else {
    logdiffexp(
      log1p(alpha) - log(alpha) + pgamma(x, alpha + 2, rate=alpha, log.p=TRUE, lower.tail=FALSE),
      2 * log(x) + pgamma(x, alpha, rate=alpha, log.p=TRUE, lower.tail=FALSE)
    ) - log(2)
  }
  if (log) {
    return(log.pev)
  } else {
    return(exp(log.pev))
  }
}

