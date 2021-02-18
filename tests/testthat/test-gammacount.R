test_that("gammacount reduces to poisson", {
  x <- seq(0, 200)
  for (lambda in c(0.01, 0.5, 1, 1.5, 10, 20, 100)) {
    expect_equal(dgc(x, lambda, alpha=1), dpois(x, lambda))
    expect_equal(pgc(x, lambda, alpha=1), ppois(x, lambda))
  }
})

test_that("gammacount cdf is sum of pmf", {
  alpha.vals <- c(0.1, 1, 10)
  lambda.vals <- c(0.5, 5, 50)
  x <- seq(0, 100)
  for (alpha in alpha.vals) {
    for (lambda in lambda.vals) {
      pmf <- dgc(x, lambda, alpha)
      cdf <- pgc(x, lambda, alpha)
      expect_equal(cumsum(pmf), cdf)
    }
  }
})

test_that("log equals the log", {
  x <- seq(0, 8)
  lambda <- rep(c(0.5, 5, 10), 3)
  alpha <- rep(c(0.1, 1, 10), 3)
  expect_equal(dgc(x, lambda, alpha, log=TRUE), log(dgc(x, lambda, alpha)))
  expect_equal(pgc(x, lambda, alpha, log.p=TRUE), log(pgc(x, lambda, alpha)))
})

test_that("upper tail equals one minus", {
  x <- seq(0, 8)
  lambda <- rep(c(0.5, 5, 10), 3)
  alpha <- rep(c(0.1, 1, 10), 3)
  expect_equal(pgc(x, lambda, alpha, lower.tail=FALSE), 1 - pgc(x, lambda, alpha))
})
