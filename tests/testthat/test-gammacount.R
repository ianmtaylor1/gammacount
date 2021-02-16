test_that("dgc reduces to dpois", {
  x <- seq(0, 200)
  for (lambda in c(0.01, 0.5, 1, 1.5, 10, 20, 100)) {
    expect_equal(dgc(x, lambda, alpha=1), dpois(x, lambda))
  }
})

test_that("pgc reduces to ppois", {
  x <- seq(0, 200)
  for (lambda in c(0.01, 0.5, 1, 1.5, 10, 20, 100)) {
    expect_equal(pgc(x, lambda, alpha=1), ppois(x, lambda))
  }
})

