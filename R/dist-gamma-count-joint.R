# This file contains functions for the joint distribution of a gamma-count with
# random start time (GCRST) random variable and the first arrival time of the
# observed count


# Distribution (density) function
dgcrst.joint <- function(x, tau, lambda, alpha, log=FALSE) {
  # Fill arguments to same length
  n <- max(length(x), length(tau), length(lambda), length(alpha))
  x <- recycle(x, n)
  tau <- recycle(tau, n)
  lambda <- recycle(lambda, n)
  alpha <- recycle(alpha, n)

  # Define the support for zero and nonzero values of x
  support.zero <- ((x == 0) & (tau > lambda))
  support.nz <- ((x > 0) & (tau <= lambda))

  # Fill in the log probability on the support. (everything else has probability
  # zero)
  lp <- rep(-Inf, n)
  lp[support.zero] <- dft(tau[support.zero], alpha[support.zero], log=TRUE)
  lp[support.nz] <- (
    dft(tau[support.nz], alpha[support.nz], log=TRUE)
    + dgc(x[support.nz] - 1, lambda[support.nz] - tau[support.nz], alpha[support.nz], log=TRUE)
  )

  if (log) {
    return(lp)
  } else {
    return(exp(lp))
  }
}


# Random number generator
rgcrst.joint <- function(n, lambda, alpha) {
  # Fill lambda to length n, in the style that other r___ function do.
  lambda.fill <- recycle(lambda, n)
  alpha.fill <- recycle(alpha, n)
  # Generate first arrival times randomly
  tau <- rft(n, alpha.fill)
  # All returned values start out at zero
  x <- rep(0, n)
  # For tau's less than lambda, generate gamma-count rv
  x[tau <= lambda.fill] <- 1 + rgc(sum(tau <= lambda.fill),
                                   lambda.fill[tau <= lambda.fill] - tau[tau <= lambda.fill],
                                   alpha.fill[tau <= lambda.fill])
  # Return both tau and x as columns in a matrix
  val <- cbind(tau, x)
  colnames(val) <- c("tau", "x")
  return(val)
}
