context("testing utils functions")

test_that("inv.logit should return values between 0 and 1", {
  expect_lt(inv.logit(10), 1)
  expect_gt(inv.logit(-10), 0)
})

test_that("partial R2 should return values between 0 and 1", {
  lm.x = lm(y_re_intercept ~ x1 + x2, data = d)
  lm.y = lm(y_re_intercept ~ x1, data = d)
  expect_lt(partialR2(lm.x, lm.y), 1)
  expect_equal(partialR2(lm.x, lm.x), 0)
  expect_gt(partialR2(lm.x, lm.y), 0)
})

test_that("binaryPGLMM within rr2 should have same results as binaryPGLMM from ape", {
  z.f1 <- ape::binaryPGLMM(y_phy_binary ~ x1, data=d, phy=phy)
  z.x1 <- ape::binaryPGLMM(y_phy_binary ~ 1, data=d, phy=phy)
  z.f <- rr2::binaryPGLMM(y_phy_binary ~ x1, data=d, phy=phy)
  z.x<- rr2::binaryPGLMM(y_phy_binary ~ 1, data=d, phy=phy)
  expect_identical(z.f1$s2, z.f$s2)
  expect_identical(z.f1$B, z.f$B)
})