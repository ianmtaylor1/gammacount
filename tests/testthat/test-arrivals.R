test_that("reduces to exponential", {
  x <- c(0.01, 0.05, seq(0.1, 10, by=.1))
  for (beta in c(0.1, 1, 10)) {
    expect_equal(dglel(x, alpha=1, beta), dexp(x, rate=beta))
    expect_equal(pglel(x, alpha=1, beta), pexp(x, rate=beta))
  }
})

test_that("reduces to gamma", {
  x <- c(0.01, 0.05, seq(0.1, 10, by=.1))
  for (beta in c(0.1, 1, 10)) {
    for (num in c(2, 5, 10, 20)) {
      expect_equal(dsga(x, num, alpha=1, beta), dgamma(x, shape=num, rate=beta))
      expect_equal(psga(x, num, alpha=1, beta), pgamma(x, shape=num, rate=beta))
    }
  }
})

test_that("pdf is derivative of cdf", {
  tol <- .Machine$double.eps ^ 0.25
  x <- c(0.01, 0.05, seq(0.1, 10, by=.1))
  alpha <- recycle(c(0.1, 1, 10, 0.1, 1, 10, 0.1, 1, 10), length(x))
  beta <- recycle(c(0.1, 0.1, 0.1, 1, 1, 1, 10, 10, 10), length(x))
  num <- recycle(c(1,2,3,4), length(x))
  eps <- .Machine$double.eps * 100

  # Take the derivative of the *log* cdf
  deriv <- numericDeriv(quote(psga(x, num, alpha, beta, log.p=TRUE)), "x")

  # What should the derivative of the *log* cdf be? d/dx log(f(x)) = f'(x)/f(x)
  expected <- exp(
    dsga(x, num, alpha, beta, log=TRUE)
    - psga(x, num, alpha, beta, log.p=TRUE)
  )

  expect_equal(diag(attr(deriv,"gradient")), expected, tolerance=tol)
})

test_that("cdf and quantile functions are inverses", {
  eps <- sqrt(.Machine$double.eps) / 2
  # p -> q
  x <- seq(0, 10, by=.1)
  alpha <- recycle(c(0.1, 1, 10), length(x))
  num <- recycle(c(1,2,3,4), length(x))
  p <- psga(x, num, alpha, log.p=TRUE, lower.tail=FALSE)
  expect_equal(qsga(p, num, alpha, log.p=TRUE, lower.tail=FALSE, tol=eps), x)
  # q -> p
  p <- seq(0.01, 0.99, by=.01)
  alpha <- recycle(c(0.1, 1, 10), length(p))
  num <- recycle(c(1,2,3,4), length(p))
  x <- qsga(p, num, alpha, tol=eps)
  expect_equal(psga(x, num, alpha), p)
})

test_that("log equals the log", {
  # density and distribution functions
  x <- seq(0, 10, by=.1)
  alpha <- recycle(c(0.1, 1, 10), length(x))
  expect_equal(dglel(x, alpha, log=TRUE), log(dglel(x, alpha)))
  expect_equal(pglel(x, alpha, log.p=TRUE), log(pglel(x, alpha)))
  # quantile functions
  p <- seq(0.01, 0.99, by=.01)
  alpha <- recycle(c(0.1, 1, 10), length(p))
  expect_equal(qglel(p, alpha), qglel(log(p), alpha, log.p=TRUE))
})

test_that("upper tail equals one minus", {
  # distribution function
  x <- seq(0, 10, by=.1)
  alpha <- recycle(c(0.1, 1, 10), length(x))
  expect_equal(pglel(x, alpha, lower.tail=FALSE), 1 - pglel(x, alpha))
  # quantile function
  p <- seq(0.01, 0.99, by=.01)
  alpha <- recycle(c(0.1, 1, 10), length(p))
  expect_equal(qglel(p, alpha), qglel(1 - p, alpha, lower.tail=FALSE))
})

test_that("output is correct length", {
  lengths <- c(3, 6, 12)
  for (alpha.length in lengths) {
    alpha <- seq_len(alpha.length)
    for (beta.length in lengths) {
      beta <- seq_len(beta.length)
      for (num.length in lengths) {
        num <- seq_len(num.length)
        for (n in lengths) {
          x <- seq_len(n)
          expect_length(rsga(n, num, alpha, beta), n)
          expected.length <- max(n, alpha.length, beta.length, num.length)
          expect_length(dsga(x, num, alpha, beta), expected.length)
          expect_length(psga(x, num, alpha, beta), expected.length)
          expect_length(qsga(x/(max(x)+1), num, alpha, beta), expected.length)
        }
      }
    }
  }
})

test_that("rng matches cdf (approximately)", {
  set.seed(123) # Set the seed so that no failures by chance
  alpha.vals <- c(0.01, 1, 10, 100)
  beta.vals <- c(0.1, 1, 10)
  num.vals <- c(1, 2, 3)
  testsig <- 0.005 # Small enough to avoid false positives
  n <- 100000 # Large enough for high power
  for (alpha in alpha.vals) {
    for (beta in beta.vals) {
      for (num in num.vals) {
        x <- rsga(n, num, alpha, beta)
        testresult <- ks.test(x, psga, num=num, alpha=alpha, beta=beta, exact=FALSE)
        expect_gt(testresult$p.value, testsig)
      }
    }
  }
})
