context("testing utils functions")

library(testthat)

test_that("inv.logit should return values between 0 and 1", {
  expect_lt(inv.logit(10), 1)
  expect_gt(inv.logit(-10), 0)
})

test_that("partial R2 should return values between 0 and 1", {
  p1 <- 10; nsample <- 10
  n <- p1 * nsample
  d <- data.frame(x1=0, x2=0, y=0, u1=rep(1:p1, each=nsample), u2=rep(1:p1, times=nsample))
  d$u1 <- as.factor(d$u1)
  d$u2 <- as.factor(d$u2)

  b1 <- 1; b2 <- -1; sd1 <- 1.5
  d$x1 <- rnorm(n=n); d$x2 <- rnorm(n=n)
  d$y <- b1 * d$x1 + b2 * d$x2 + rep(rnorm(n=p1, sd=sd1), each=nsample) + rep(rnorm(n=p1, sd=sd1), times=nsample) + rnorm(n=n)
  lm.x = lm(y ~ x1 + x2, data = d)
  lm.y = lm(y ~ x1, data = d)
  expect_lt(partialR2(lm.x, lm.y), 1)
  expect_equal(partialR2(lm.x, lm.x), 0)
  expect_gt(partialR2(lm.x, lm.y), 0)
})
