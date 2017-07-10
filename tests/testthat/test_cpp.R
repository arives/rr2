context("testing R2.ce functions with loops written in c++")

test_that("PGLS: c++ version give the same results with the R version", {
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

  expect_equal(rr2::R2.ce(mod, mod.r, phy = phy, cpp = T), rr2::R2.ce(mod, mod.r, phy = phy, cpp = F))
  
  # try again
  p1 <- 11; nsample <- 11; n <- p1 * nsample
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
  
  expect_equal(rr2::R2.ce(mod, mod.r, phy = phy, cpp = T), rr2::R2.ce(mod, mod.r, phy = phy, cpp = F))
})