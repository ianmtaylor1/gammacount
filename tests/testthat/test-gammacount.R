test_that("reduces to poisson", {
  x <- seq(0, 200)
  for (lambda in c(0.1, 0.5, 1, 1.5, 10)) {
    for (beta in c(0.1, 0.5, 1, 1.5, 10)) {
      expect_equal(dgc(x, lambda, alpha=1, beta), dpois(x, lambda * beta))
      expect_equal(pgc(x, lambda, alpha=1, beta), ppois(x, lambda * beta))
    }
  }
})

test_that("cdf is sum of pmf", {
  alpha.vals <- c(0.1, 1, 10)
  lambda.vals <- c(0.5, 5, 10)
  beta.vals <- c(0.5, 5, 10)
  x <- seq(0, 100)
  for (alpha in alpha.vals) {
    for (lambda in lambda.vals) {
      for (beta in beta.vals) {
        pmf <- dgc(x, lambda, alpha, beta)
        cdf <- pgc(x, lambda, alpha, beta)
        expect_equal(cumsum(pmf), cdf)
      }
    }
  }
})

test_that("cdf and quantile functions are inverses", {
  # Test p -> q
  # Uses lower.tail=F and log.p=T to handle large values of x gracefully
  x <- seq(0, 100)
  lambda <- recycle(c(0.5, 5, 10, 0.5, 5, 10, 0.5, 5, 10), length(x))
  alpha <- recycle(c(0.1, 0.1, 0.1, 1, 1, 1, 10, 10, 10), length(x))
  p <- pgc(x, lambda, alpha, lower.tail=FALSE, log.p=TRUE)
  expect_equal(qgc(p, lambda, alpha, lower.tail=FALSE, log.p=TRUE), x)
  # Test q -> p
  p <- seq(0.01, 0.99, by=.01)
  lambda <- recycle(c(0.5, 5, 10, 0.5, 5, 10, 0.5, 5, 10), length(p))
  alpha <- recycle(c(0.1, 0.1, 0.1, 1, 1, 1, 10, 10, 10), length(p))
  q <- qgc(p, lambda, alpha)
  expect_true(all(pgc(q, lambda, alpha) >= p))
  expect_true(all(pgc(q-1, lambda, alpha) < p))
})

test_that("log equals the log", {
  # density and distribution functions
  x <- seq(0, 8)
  lambda <- recycle(c(0.5, 5, 10, 0.5, 5, 10, 0.5, 5, 10), length(x))
  alpha <- recycle(c(0.1, 0.1, 0.1, 1, 1, 1, 10, 10, 10), length(x))
  expect_equal(dgc(x, lambda, alpha, log=TRUE), log(dgc(x, lambda, alpha)))
  expect_equal(pgc(x, lambda, alpha, log.p=TRUE), log(pgc(x, lambda, alpha)))
  # quantile functions
  p <- seq(0.01, 0.99, by=.01)
  lambda <- recycle(c(0.5, 5, 10, 0.5, 5, 10, 0.5, 5, 10), length(p))
  alpha <- recycle(c(0.1, 0.1, 0.1, 1, 1, 1, 10, 10, 10), length(p))
  expect_equal(qgc(p, lambda, alpha), qgc(log(p), lambda, alpha, log.p=TRUE))
})

test_that("upper tail equals one minus", {
  # distribution function
  x <- seq(0, 8)
  lambda <- recycle(c(0.5, 5, 10, 0.5, 5, 10, 0.5, 5, 10), length(x))
  alpha <- recycle(c(0.1, 0.1, 0.1, 1, 1, 1, 10, 10, 10), length(x))
  expect_equal(pgc(x, lambda, alpha, lower.tail=FALSE), 1 - pgc(x, lambda, alpha))
  # quantile function
  p <- seq(0.01, 0.99, by=.01)
  lambda <- recycle(c(0.5, 5, 10, 0.5, 5, 10, 0.5, 5, 10), length(p))
  alpha <- recycle(c(0.1, 0.1, 0.1, 1, 1, 1, 10, 10, 10), length(p))
  expect_equal(qgc(p, lambda, alpha), qgc(1 - p, lambda, alpha, lower.tail=FALSE))
})

test_that("output is correct length", {
  lengths <- c(3, 6, 12)
  for (alpha.length in lengths) {
    alpha <- seq_len(alpha.length)
    for (lambda.length in lengths) {
      lambda <- seq_len(lambda.length)
      for (beta.length in lengths) {
        beta <- seq_len(beta.length)
        for (n in lengths) {
          x <- seq_len(n)
          expect_length(rgc(n, lambda, alpha, beta), n)
          expected.length <- max(n, lambda.length, alpha.length, beta.length)
          expect_length(dgc(x, lambda, alpha, beta), expected.length)
          expect_length(pgc(x, lambda, alpha, beta), expected.length)
          expect_length(qgc(x/(max(x)+1), lambda, alpha, beta), expected.length)
        }
      }
    }
  }
})
