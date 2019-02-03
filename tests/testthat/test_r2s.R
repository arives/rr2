context("testing R2 functions")

# data 
set.seed(123456)
p1 <- 10; nsample <- 20; n <- p1 * nsample
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

# data for communityPGLMM
# library(dplyr)
# comm = phyr::comm_a
# comm$site = row.names(comm)
# dat = tidyr::gather(comm, key = "sp", value = "freq", -site) %>%
#   dplyr::left_join(phyr::envi, by = "site") %>%
#   dplyr::left_join(phyr::traits, by = "sp")
# dat$pa = as.numeric(dat$freq > 0)

test_that("when missing mod.r, the functions will automatically create one", {
  testthat::skip_on_cran() # don't run on CRAN since these are time consuming
  # LM
  z.f <- lm(y_re_intercept ~ x1 + x2, data = d)
  z.0 <- lm(y_re_intercept ~ 1, data = d)
  expect_equal(R2(z.f, z.0), R2(z.f))
  
  # GLM
  z.f.glm <- glm(y_binary ~ x1 + x2, data = d, family = "binomial")
  z.v.glm <- glm(y_binary ~ 1, data = d, family = "binomial")
  expect_equal(R2(mod = z.f.glm, mod.r = z.v.glm),
               R2(mod = z.f.glm))
  
  # LMM
  z.f <- lme4::lmer(y_re_intercept ~ x1 + x2 + (1 | u1) + (1 | u2), data = d, REML = F)
  z.0 <- lm(y_re_intercept ~ 1, data = d)
  expect_equal(R2(z.f, z.0), R2(z.f))
  
  # GLMM
  z.f.glmm <- lme4::glmer(y_binary ~ x1 + (1 | u1), data = d, family = "binomial")
  z.v.glmm <- glm(y_binary ~ 1, data = d, family = "binomial")
  expect_equal(R2(mod = z.f.glmm, mod.r = z.v.glmm),
               R2(mod = z.f.glmm))
})

test_that("when missing mod.r, the functions will automatically create one, long-run analyses", {
  testthat::skip_on_cran() # don't run on CRAN since these are time consuming
  
  # PGLS
  z.x.pgls <- phylolm::phylolm(y_pgls ~ 1, phy = phy, data = d, model = "lambda")
  lam.x.pgls <- round(z.x.pgls$optpar, digits = 4)
  z.f.pgls <- phylolm::phylolm(y_pgls ~ x_trait, phy = phy, data = d, model = "lambda")
  z.v.pgls <- lm(y_pgls ~ 1, data = d)
  expect_equal(R2(mod = z.f.pgls, mod.r = z.v.pgls, phy = phy),
               R2(mod = z.f.pgls, phy = phy))
  
  for(md in c("OUrandomRoot", "OUfixedRoot", "BM", 
              "kappa", "delta", "EB")){
    z.f.pgls2 <- phylolm::phylolm(y_pgls ~ x_trait, phy = phy, data = d, model = md)
    z.v.pgls2 <- lm(y_pgls ~ 1, data = d)
    expect_equal(R2(mod = z.f.pgls2, mod.r = z.v.pgls2, phy = phy),
                 R2(mod = z.f.pgls2, phy = phy))
  }
  
  # GLS
  z.f.gls <- nlme::gls(y_pgls ~ x_trait, data = d, correlation = ape::corPagel(1, phy), method = "ML")
  expect_equal(R2(mod = z.f.gls, mod.r = z.v.pgls),
               R2(mod = z.f.gls))
  
  # binaryPGLMM
  z.f.plog <- rr2::binaryPGLMM(y_phy_binary ~ x1, data = d, phy = phy)
  z.v.plog <- glm(y_phy_binary ~ 1, data = d, family = "binomial")
  # R2.lik can't be used with binaryPGLMM because it is not a ML method
  expect_message(t1 <- R2(mod = z.f.plog, mod.r = z.v.plog), 
                 "Models of class binaryPGLMM do not have R2.lik method")
  expect_message(t2 <- R2(mod = z.f.plog, mod.r = z.v.plog), 
                 "Models of class binaryPGLMM do not have R2.lik method")
  expect_equal(t1, t2)
  
  z.f.plog2 <- phylolm::phyloglm(y_phy_binary ~ x1, data = d, start.alpha = 1, phy = phy)
  z.v.plog2 <- glm(y_phy_binary ~ 1, data = d, family = "binomial")
  expect_message(t11 <- R2(z.f.plog2, z.v.plog2),
                 "Models of class phyloglm only have R2.lik method")
  expect_message(t12 <- R2(z.f.plog2),
                 "Models of class phyloglm only have R2.lik method")
  expect_equal(t11, t12)
  
  # # communityPGLMM
  # mod <- phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),
  #                                           dat, tree = phyr::phylotree, REML = F)
  # mod.r <- lm(freq ~ 1, dat)
  # expect_message(R2(mod, mod.r))
  # expect_equal(R2(mod, mod.r), R2(mod))
})


test_that("when lmer models were fitted with REML = T, 
          R2.lr (but not R2.ls and R2.ce) will change it to FALSE", {
  z.f2 <- lme4::lmer(y_re_intercept ~ x1 + x2 + (1 | u1) + (1 | u2), data = d, REML = T)
  expect_warning(R2(z.f2), "mod updated with REML = F")
})
