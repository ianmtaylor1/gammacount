# This file contains functions for the joint distribution of a gamma-count with
# random start time (GCRST) random variable and the first arrival time of the
# observed count


#' The Joint Gamma-Count RST and First Arrival Time Distribution
#'
#' Density, and random number generation for the joint gamma-count distribution
#' with random start time and first arrival time, with scale/mean parameter
#' lambda and dispersion parameter alpha
#'
#' The gamma-count distribution with random start time models the count of
#' events arriving in the interval (L,L+lambda] where L is a large randomly
#' chosen positive real number and the times between successive events follow a
#' gamma(alpha, alpha) distribution. The first arrival time is the time between
#' L and the first event which arrives after L, and may be greater than
#' L+lambda. The joint distribution of both (arising from the same process) is
#' given by these functions.
#'
#' @param x vector of (non-negative integer) quantiles
#' @param tau vector of (positive) first arrival time values
#' @param n number of random values to return
#' @param lambda vector of (positive) scale/mean parameters
#' @param alpha vector of (positive) dispersion parameters
#' @param log logical; if TRUE, probabilities p are given as log(p)
#'
#' @return dgcrst.joint gives the (log) joint density, and rgcrst.joint
#'   generates random counts and arrival times.
#'
#'   Invalid lambda or alpha will result in return value NaN, with a warning.
#'
#'   The length of the result is determined by n for rgcrst.joint, and the
#'   length of the longest numeric argument for the other functions. The
#'   numerical arguments other than n are recycled to the length of the result.
#'   Only the first elements of the logical arguments are used.
#'
#'   In dgcrst and pgcrst, if diagnostics=TRUE, the returned value is a list with
#'   components value, abs.error, subdivisions, message. All components are
#'   vectors with the same length. 'value' is the probability or density
#'   requested, and all other values come directly from any calls
#'   to integrate(). NA is put in the diagnostics for any value for which
#'   integrate() did not need to be called.
#'
#' @seealso [dgc()] for the standard gamma-count distribution, and [dft()] for
#'   the first arrival time distribution after the random start time, and
#'   [dgcrst()] for the count distribution with random start time.
#'
#' @name GammaCountRSTJoint
NULL


# Distribution (density) function
#' @rdname GammaCountRSTJoint
#' @export
dgcrst.joint <- function(x, tau, lambda, alpha=1, log=FALSE) {
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
#' @rdname GammaCountRSTJoint
#' @export
rgcrst.joint <- function(n, lambda, alpha=1) {
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
  val <- cbind(x, tau)
  colnames(val) <- c("x", "tau")
  return(val)
}
