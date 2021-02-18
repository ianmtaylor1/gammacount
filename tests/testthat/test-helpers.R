
#### Tests of logdiffexp


test_that("valid finite logdiffexp inputs", {
  # Test various correct (i.e. a > b) finite inputs for numerical accuracy
  base <- c(-100, 0, 100)
  inc <- c(0, log(2)/2, log(2)+2, 100)
  a <- c(outer(base, inc, "+"))
  b <- c(outer(base, rep(0, length(inc)), "+"))
  expect_equal(logdiffexp(a, b), log(exp(a) - exp(b)))
})

test_that("invalid finite logdiffexp inputs", {
  # Test that a warning is issued when b > a
  expect_warning(logdiffexp(0, 0.001), "NaNs produced")
  expect_warning(logdiffexp(0, 100), "NaNs produced")
})

test_that("infinite logdiffexp inputs", {
  # Test valid (a > b) inputs involving at least one infinite value
  a.valid <- c(-Inf, 0,    Inf,  Inf)
  b.valid <- c(-Inf, -Inf, -Inf, 0)
  expect_equal(logdiffexp(a.valid, b.valid), log(exp(a.valid) - exp(b.valid)))

  # Test infinite inputs that should result in a warning and NaNs
  a.warn <- c(Inf, -Inf, -Inf, 0)
  b.warn <- c(Inf, 0, Inf, Inf)
  for (i in 1:length(a.warn)) {
    expect_warning(logdiffexp(a.warn[i], b.warn[i]), "NaNs produced")
  }
})

test_that("missing logdiffexp inputs", {
  a <- c(NaN, -Inf, 0, Inf, NaN, NaN, NaN)
  b <- c(NaN, NaN, NaN, NaN, -Inf, 0, Inf)
  expect_equal(logdiffexp(a, b), rep(NaN, length(a)))
})

#### Tests of logsumexp

test_that("valid logsumexp inputs", {
  a <- c(Inf, Inf, -Inf,
         Inf, -Inf,
         -100, 0, 100, 100)
  b <- c(Inf, -Inf, -Inf,
         0, 0,
         -100, 0, 100, -100)
  expect_equal(logsumexp(a, b), log(exp(a) + exp(b)))
  expect_equal(logsumexp(b, a), log(exp(b) + exp(a)))
})

test_that("missing logsumexp inputs", {
  a <- c(NaN, -Inf, 0, Inf, NaN, NaN, NaN)
  b <- c(NaN, NaN, NaN, NaN, -Inf, 0, Inf)
  expect_equal(logsumexp(a, b), rep(NaN, length(a)))
})

#### Tests of recycle

test_that("recycle to correct length", {
  v <- c(1,2,3,4,5)
  for (n in c(1, 3, 5, 7, 10, 102)) {
    expect_length(recycle(v, n), n)
  }
})
