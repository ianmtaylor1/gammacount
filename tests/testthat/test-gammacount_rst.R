test_that("reduces to poisson", {
  tol <- .Machine$double.eps ^ 0.25
  x <- seq(0, 200)
  for (lambda in c(0.01, 0.5, 1, 1.5, 10, 20, 100)) {
    expect_equal(dgcrst(x, lambda, alpha=1), dpois(x, lambda), tolerance=tol)
    expect_equal(pgcrst(x, lambda, alpha=1), ppois(x, lambda), tolerance=tol)
  }
})

test_that("cdf is sum of pmf", {
  tol <- .Machine$double.eps ^ 0.25
  alpha.vals <- c(0.1, 1, 10)
  lambda.vals <- c(0.5, 5, 50)
  x <- seq(0, 100)
  for (alpha in alpha.vals) {
    for (lambda in lambda.vals) {
      pmf <- dgcrst(x, lambda, alpha)
      cdf <- pgcrst(x, lambda, alpha)
      expect_equal(cdf, cumsum(pmf), tolerance=tol)
    }
  }
})

test_that("cdf and quantile functions are inverses", {
  alpha.vals <- c(0.1, 1, 10)
  lambda.vals <- c(0.5, 5, 10)
  for (alpha in alpha.vals) {
    for (lambda in lambda.vals) {
      x <- seq(0, ceiling(1.9 * lambda))
      # Test p -> q
      p <- pgcrst(x, lambda, alpha)
      expect_equal(qgcrst(p, lambda, alpha), x)
      # Test q -> p
      p <- seq(0.01, 0.99, by=.01)
      q <- qgcrst(p, lambda, alpha)
      expect_true(all(pgcrst(q, lambda, alpha) >= p))
      expect_true(all(pgcrst(q-1, lambda, alpha) < p))
    }
  }
})

test_that("log equals the log", {
  # density and distribution functions
  x <- seq(0, 8)
  lambda <- recycle(c(0.5, 5, 10, 0.5, 5, 10, 0.5, 5, 10), length(x))
  alpha <- recycle(c(0.1, 0.1, 0.1, 1, 1, 1, 10, 10, 10), length(x))
  expect_equal(dgcrst(x, lambda, alpha, log=TRUE), log(dgcrst(x, lambda, alpha)))
  expect_equal(pgcrst(x, lambda, alpha, log.p=TRUE), log(pgcrst(x, lambda, alpha)))
  # quantile functions
  p <- seq(0.01, 0.99, by=.01)
  lambda <- recycle(c(0.5, 5, 10, 0.5, 5, 10, 0.5, 5, 10), length(p))
  alpha <- recycle(c(0.1, 0.1, 0.1, 1, 1, 1, 10, 10, 10), length(p))
  expect_equal(qgcrst(p, lambda, alpha), qgcrst(log(p), lambda, alpha, log.p=TRUE))
})

test_that("upper tail equals one minus", {
  # distribution function
  x <- seq(0, 8)
  lambda <- recycle(c(0.5, 5, 10, 0.5, 5, 10, 0.5, 5, 10), length(x))
  alpha <- recycle(c(0.1, 0.1, 0.1, 1, 1, 1, 10, 10, 10), length(x))
  expect_equal(pgcrst(x, lambda, alpha, lower.tail=FALSE), 1 - pgcrst(x, lambda, alpha))
  # quantile function
  p <- seq(0.01, 0.99, by=.01)
  lambda <- recycle(c(0.5, 5, 10, 0.5, 5, 10, 0.5, 5, 10), length(p))
  alpha <- recycle(c(0.1, 0.1, 0.1, 1, 1, 1, 10, 10, 10), length(p))
  expect_equal(qgcrst(p, lambda, alpha), qgcrst(1 - p, lambda, alpha, lower.tail=FALSE))
})

test_that("output is correct length", {
  lengths <- c(3, 6, 12)
  for (alpha.length in lengths) {
    alpha <- seq_len(alpha.length)
    for (lambda.length in lengths) {
      lambda <- seq_len(lambda.length)
      for (n in lengths) {
        x <- seq_len(n)
        expect_length(rgcrst(n, lambda, alpha), n)
        expected.length <- max(n, lambda.length, alpha.length)
        expect_length(dgcrst(x, lambda, alpha), expected.length)
        expect_length(pgcrst(x, lambda, alpha), expected.length)
        expect_length(qgcrst(x/(max(x)+1), lambda, alpha), expected.length)
      }
    }
  }
})
