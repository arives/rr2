context("testing R2 functions")

# data 
# set.seed(123)
p1 <- 10; nsample <- 10; n <- p1 * nsample
d <- data.frame(x1 = rnorm(n = n), 
                x2 = rnorm(n = n), 
                u1 = rep(1:p1, each = nsample), 
                u2 = rep(1:p1, times = nsample))
d$u1 <- as.factor(d$u1); d$u2 <- as.factor(d$u2)

# LMM: y with random intercept
b1 <- 1; b2 <- -1; sd1 <- 1.5
d$y_re_intercept <- b1 * d$x1 + b2 * d$x2 + 
  rep(rnorm(n = p1, sd = sd1), each = nsample) +  # random intercept u1
  rep(rnorm(n = p1, sd = sd1), times = nsample) + # random intercept u2
  rnorm(n = n)

# LMM: y with random slope
b1 <- 0; sd1 <- 1; sd.x1 <- 2
d$y_re_slope <- b1 * d$x1 + 
  rep(rnorm(n = p1, sd = sd1), each = nsample) + # random intercept u1
  d$x1 * rep(rnorm(n = p1, sd = sd.x1), times = nsample) + # random slope u1
  rnorm(n = n)

# GLMM
b1 <- 1; sd1 <- 1.5
prob <- rr2::inv.logit(b1 * d$x1 + rep(rnorm(n = p1, sd = sd1), each = nsample)) # random intercept u1
d$y_binary <- rbinom(n = n, size = 1, prob = prob)

# y for PGLS
b1 <- 1.5; signal <- 0.7
phy <- ape::compute.brlen(ape::rtree(n = n), method = "Grafen", power = 1)
phy.x <- ape::compute.brlen(phy, method = "Grafen", power = .0001)
x_trait <- ape::rTraitCont(phy.x, model = "BM", sigma = 1)
e <- signal^0.5 * ape::rTraitCont(phy, model = "BM", sigma = 1) + (1-signal)^0.5 * rnorm(n=n)
d$x_trait <- x_trait[match(names(e), names(x_trait))]
d$y_pgls <- b1 * x_trait + e
rownames(d) <- phy$tip.label	

# y for Phylogenetic logistic regression
b1 <- 1.5; signal <- 2
e <- signal * ape::rTraitCont(phy, model = "BM", sigma = 1)
e <- e[match(phy$tip.label, names(e))]
d$y_phy_binary <- rbinom(n = n, size = 1, prob = rr2::inv.logit(b1 * d$x1 + e))

test_that("R2 functions do not work with lm", {
    z.lm = lm(y_re_intercept ~ x1 + x2, data = d)
    expect_error(R2.lr(z.lm), "mod must be class one of classes lmerMod, glmerMod, phylolm, phyloglm.")
})

test_that("when missing mod.r, the functions will automatically creat one", {
    z.f <- lme4::lmer(y_re_intercept ~ x1 + x2 + (1 | u1) + (1 | u2), data = d, REML = F)
    z.0 <- lm(y_re_intercept ~ 1, data = d)
    
    expect_equal(R2.ls(z.f, z.0), R2.ls(z.f))
    expect_equal(R2.ce(z.f, z.0), R2.ce(z.f))
    expect_equal(R2.lr(z.f, z.0), R2.lr(z.f))
})

test_that("when lmer models were fitted with REML = T, R2.lr (but not R2.ls and R2.ce) will change it to FALSE", {
    z.f2 <- lme4::lmer(y_re_intercept ~ x1 + x2 + (1 | u1) + (1 | u2), data = d, REML = T)
    expect_warning(R2.lr(z.f2), "mod updated with REML=F")
})
