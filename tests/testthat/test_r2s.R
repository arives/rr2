context("testing R2 functions")

set.seed(123)
p1 <- 10; nsample <- 10
n <- p1 * nsample
d <- data.frame(x1=0, x2=0, y=0, u1=rep(1:p1, each=nsample), u2=rep(1:p1, times=nsample))
d$u1 <- as.factor(d$u1)
d$u2 <- as.factor(d$u2)
b1 <- 1; b2 <- -1; sd1 <- 1.5
d$x1 <- rnorm(n=n); d$x2 <- rnorm(n=n)
d$y <- b1 * d$x1 + b2 * d$x2 + rep(rnorm(n=p1, sd=sd1), each=nsample) + rep(rnorm(n=p1, sd=sd1), times=nsample) + rnorm(n=n)

test_that("R2 functions do not work with lm", {
  z.lm = lm(y ~ x1 + x2, data = d)
  expect_error(R2.lr(z.lm), "mod must be class one of classes lmerMod, glmerMod, phylolm, phyloglm.")
})

test_that("when missing mod.r, the functions will automatically creat one", {
  z.f <- lme4::lmer(y ~ x1 + x2 + (1 | u1) + (1 | u2), data=d, REML = F)
  z.0 <- lm(y ~ 1, data=d)

  expect_equal(R2.ls(z.f, z.0), R2.ls(z.f))
  expect_equal(R2.ce(z.f, z.0), R2.ce(z.f))
  expect_equal(R2.lr(z.f, z.0), R2.lr(z.f))
})

test_that("when lmer models were fitted with REML=T, change it to F", {
  z.f2 <- lme4::lmer(y ~ x1 + x2 + (1 | u1) + (1 | u2), data=d, REML = T)
  expect_warning(R2.lr(z.f2, z.0), "mod updated with REML=F")
})




