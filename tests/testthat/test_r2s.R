context("testing R2 functions")

test_that("R2 functions do not work with lm", {
  z.lm = lm(y_re_intercept ~ x1 + x2, data = d)
  expect_error(R2.lr(z.lm), "mod must be class one of classes lmerMod, glmerMod, phylolm, phyloglm.")
})

test_that("when missing mod.r, the functions will automatically creat one", {
  z.f <- lme4::lmer(y_re_intercept ~ x1 + x2 + (1 | u1) + (1 | u2), data=d, REML = F)
  z.0 <- lm(y_re_intercept ~ 1, data=d)

  expect_equal(R2.ls(z.f, z.0), R2.ls(z.f))
  expect_equal(R2.ce(z.f, z.0), R2.ce(z.f))
  expect_equal(R2.lr(z.f, z.0), R2.lr(z.f))
})

test_that("when lmer models were fitted with REML = T, R2.lr (but not R2.ls and R2.ce) will change it to FALSE", {
  z.f2 <- lme4::lmer(y_re_intercept ~ x1 + x2 + (1 | u1) + (1 | u2), data = d, REML = T)
  expect_warning(R2.lr(z.f2), "mod updated with REML=F")
})
