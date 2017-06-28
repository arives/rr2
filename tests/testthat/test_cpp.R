context("testing R2.ce functions with loops written in c++")

test_that("lmer: c++ version give the same results with the R version", {
  set.seed(123)
  p1 <- 10; nsample <- 10
  n <- p1 * nsample
  d <- data.frame(x1=0, x2=0, y=0, u1=rep(1:p1, each=nsample), u2=rep(1:p1, times=nsample))
  d$u1 <- as.factor(d$u1)
  d$u2 <- as.factor(d$u2)
  b1 <- 1; b2 <- -1; sd1 <- 1.5
  d$x1 <- rnorm(n=n); d$x2 <- rnorm(n=n)
  d$y <- b1 * d$x1 + b2 * d$x2 + rep(rnorm(n=p1, sd=sd1), each=nsample) + rep(rnorm(n=p1, sd=sd1), times=nsample) + rnorm(n=n)

  mod = lme4::lmer(y ~ x1 + x2 + (1 | u1) + (1 | u2), data=d, REML = F)
  mod.r1 = lme4::lmer(y ~ x1 + (1 | u1) + (1 | u2), data=d, REML = F)
  mod.r2 = lme4::lmer(y ~ 1 + (1 | u1) + (1 | u2), data=d, REML = F)
  mod.r3 = lm(y ~ x1, data=d)

  expect_equal(R2.ce(mod, mod.r1, cpp = T), R2.ce(mod, mod.r1, cpp = F))
  expect_equal(R2.ce(mod, mod.r2, cpp = T), R2.ce(mod, mod.r2, cpp = F))
  expect_equal(R2.ce(mod, mod.r3, cpp = T), R2.ce(mod, mod.r3, cpp = F))
})

test_that("glmer: c++ version give the same results with the R version", {
  set.seed(123)
  p1 <- 10; nsample <- 10; n <- p1 * nsample
  d <- data.frame(x=0, y=0, u=rep(1:p1, each=nsample))
  d$u <- as.factor(d$u)
  b1 <- 1; sd1 <- 1.5
  d$x <- rnorm(n=n)
  prob <- inv.logit(b1 * d$x + rep(rnorm(n=p1, sd=sd1), each=nsample))
  d$y <- rbinom(n=n, size=1, prob=prob)

  mod <- lme4::glmer(y ~ x + (1 | u), data=d, family="binomial")
  mod.r1 <- lme4::glmer(y ~ 1 + (1 | u), data=d, family="binomial")
  mod.r2 <- glm(y ~ x, data=d, family="binomial")

  expect_equal(R2.ce(mod, mod.r1, cpp = T), R2.ce(mod, mod.r1, cpp = F))
  expect_equal(R2.ce(mod, mod.r2, cpp = T), R2.ce(mod, mod.r2, cpp = F))
})

test_that("PGLS: c++ version give the same results with the R version", {
  set.seed(123)
  p1 <- 10; nsample <- 10; n <- p1 * nsample
  d <- data.frame(x=array(0, dim=n), y=0)
  b1 <- 1.5; signal <- 0.7
  phy <- ape::compute.brlen(ape::rtree(n=n), method = "Grafen", power = 1)
  phy.x <- phy
  x <- ape::rTraitCont(phy.x, model = "BM", sigma = 1)
  e <- signal^0.5 * ape::rTraitCont(phy, model = "BM", sigma = 1) + (1-signal)^0.5 * rnorm(n=n)
  d$x <- x[match(names(e), names(x))]
  d$y <- b1 * x + e
  rownames(d) <- phy$tip.label

  z.x <- phylolm::phylolm(y ~ 1, phy=phy, data=d, model="lambda")
  lam.x <- round(z.x$optpar, digits=4)
  mod <- phylolm::phylolm(y ~ x, phy=phy, data=d, model="lambda", starting.value=.98*lam.x+.01)
  mod.r <- lm(y ~ x, data=d)

  expect_equal(R2.ce(mod, mod.r, phy = phy, cpp = T), R2.ce(mod, mod.r, phy = phy, cpp = F))
})

test_that("PLOG: c++ version give the same results with the R version", {
  set.seed(123)
  p1 <- 10; nsample <- 10; n <- p1 * nsample
  b1 <- 1.5; signal <- 2
  phy <- ape::compute.brlen(ape::rtree(n=n), method = "Grafen", power = 1)
  d <- data.frame(x=array(0, dim=n), y=0)
  d$x <- rnorm(n)
  e <- signal * ape::rTraitCont(phy, model = "BM", sigma = 1)
  e <- e[match(phy$tip.label, names(e))]
  d$y <- rbinom(n=n, size=1, prob=rr2::inv.logit(b1 * d$x + e))
  rownames(d) <- phy$tip.label

  mod <- rr2::binaryPGLMM(formula = "y ~ x", data=d, phy=phy)
  mod.r1 <- rr2::binaryPGLMM(y ~ 1, data=d, phy=phy)
  mod.r2 <- glm(y ~ x, data=d, family="binomial")

  expect_equal(R2.ce(mod, mod.r1, cpp = T), R2.ce(mod, mod.r1, cpp = F))
  expect_equal(R2.ce(mod, mod.r2, cpp = T), R2.ce(mod, mod.r2, cpp = F))
})
