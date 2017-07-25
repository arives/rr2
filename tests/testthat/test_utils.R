context("testing utils functions")

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

# there is a typo in ape::binaryPGLMM, which was corrected in rr2::binaryPGLMM, so not equal now...
# test_that("binaryPGLMM within rr2 should have same results as binaryPGLMM from ape", {
#     z.f1 <- ape::binaryPGLMM(y_phy_binary ~ x1, data = d, phy = phy)
#     z.x1 <- ape::binaryPGLMM(y_phy_binary ~ 1, data = d, phy = phy)
#     z.f <- rr2::binaryPGLMM(y_phy_binary ~ x1, data = d, phy = phy)
#     z.x <- rr2::binaryPGLMM(y_phy_binary ~ 1, data = d, phy = phy)
#     expect_identical(z.f1$s2, z.f$s2)
#     expect_identical(z.f1$B, z.f$B)
# })
