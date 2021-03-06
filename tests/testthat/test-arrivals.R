test_that("reduces to exponential", {
  x <- c(0.01, 0.05, seq(0.1, 10, by=.1))
  expect_equal(dft(x, alpha=1), dexp(x, rate=1))
  expect_equal(pft(x, alpha=1), pexp(x, rate=1))
})

test_that("pdf is derivative of cdf", {
  tol <- .Machine$double.eps ^ 0.25
  x <- c(0.01, 0.05, seq(0.1, 10, by=.1))
  alpha <- recycle(c(0.1, 1, 10), length(x))
  eps <- .Machine$double.eps * 100
  deriv <- numericDeriv(quote(pft(x, alpha)), "x")
  expect_equal(diag(attr(deriv,"gradient")), dft(x, alpha), tolerance=tol)
})

test_that("cdf and quantile functions are inverses", {
  eps <- sqrt(.Machine$double.eps) / 2
  # p -> q
  x <- seq(0, 10, by=.1)
  alpha <- recycle(c(0.1, 1, 10), length(x))
  p <- pft(x, alpha, log.p=TRUE, lower.tail=FALSE)
  expect_equal(qft(p, alpha, log.p=TRUE, lower.tail=FALSE, tol=eps), x)
  # q -> p
  p <- seq(0.01, 0.99, by=.01)
  alpha <- recycle(c(0.1, 1, 10), length(p))
  x <- qft(p, alpha, tol=eps)
  expect_equal(pft(x, alpha), p)
})

test_that("log equals the log", {
  # density and distribution functions
  x <- seq(0, 10, by=.1)
  alpha <- recycle(c(0.1, 1, 10), length(x))
  expect_equal(dft(x, alpha, log=TRUE), log(dft(x, alpha)))
  expect_equal(pft(x, alpha, log.p=TRUE), log(pft(x, alpha)))
  # quantile functions
  p <- seq(0.01, 0.99, by=.01)
  alpha <- recycle(c(0.1, 1, 10), length(p))
  expect_equal(qft(p, alpha), qft(log(p), alpha, log.p=TRUE))
})

test_that("upper tail equals one minus", {
  # distribution function
  x <- seq(0, 10, by=.1)
  alpha <- recycle(c(0.1, 1, 10), length(x))
  expect_equal(pft(x, alpha, lower.tail=FALSE), 1 - pft(x, alpha))
  # quantile function
  p <- seq(0.01, 0.99, by=.01)
  alpha <- recycle(c(0.1, 1, 10), length(p))
  expect_equal(qft(p, alpha), qft(1 - p, alpha, lower.tail=FALSE))
})

test_that("output is correct length", {
  lengths <- c(3, 6, 12)
  for (alpha.length in lengths) {
    alpha <- seq_len(alpha.length)
    for (n in lengths) {
      x <- seq_len(n)
      expect_length(rft(n, alpha), n)
      expected.length <- max(n, alpha.length)
      expect_length(dft(x, alpha), expected.length)
      expect_length(pft(x, alpha), expected.length)
      expect_length(qft(x/(max(x)+1), alpha), expected.length)
    }
  }
})

test_that("rng matches cdf (approximately)", {
  set.seed(123) # Set the seed so that no failures by chance
  alpha.vals <- c(0.01, 0.1, 1, 10, 100)
  testsig <- 0.01 # Small enough to avoid false positives
  n <- 100000 # Large enough for high power
  for (alpha in alpha.vals) {
    x <- rft(n, alpha)
    testresult <- ks.test(x, pft, alpha=alpha)
    expect_gt(testresult$p.value, testsig)
  }
})

test_that("partial ev passes basic checks", {
  # Check that full calculation leads to EV.
  alpha <- c(0.01, 0.1, 1, 10, 100)
  expect_equal(ft.partial.ev(0, alpha, lower.tail=FALSE), (alpha + 1) / (2 * alpha))
  # Check that other boundaries return zero
  expect_equal(ft.partial.ev(0, alpha), rep(0, length(alpha)))
  # Lower tail and upper tail add up to expected value
  x <- seq(0.1, 10, by=0.1)
  alpha <- recycle(alpha, length(x))
  expect_equal(
    ft.partial.ev(x, alpha, lower.tail=TRUE) +
      ft.partial.ev(x, alpha, lower.tail=FALSE),
    (alpha + 1) / (2 * alpha))
})
