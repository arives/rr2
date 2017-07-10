context("testing R2.ce functions with loops written in c++")

test_that("PGLS: c++ version give the same results with the R version", {
  z.x <- phylolm::phylolm(y_pgls ~ 1, phy = phy, data = d, model = "lambda")
  lam.x <- round(z.x$optpar, digits=4)
  mod <- phylolm::phylolm(y_pgls ~ x_trait, phy=phy, data=d, model="lambda", starting.value=.98*lam.x+.01)
  mod.r <- lm(y_pgls ~ x_trait, data=d)

  expect_equal(rr2::R2.ce(mod, mod.r, phy = phy, cpp = T), rr2::R2.ce(mod, mod.r, phy = phy, cpp = F))
})