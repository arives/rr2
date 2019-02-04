<!-- README.md is generated from README.Rmd. Please edit that file -->
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01028/status.svg)](https://doi.org/10.21105/joss.01028)
[![CRAN
status](https://www.r-pkg.org/badges/version/rr2)](https://cran.r-project.org/package=rr2)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/rr2)](http://www.r-pkg.org/pkg/rr2)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/grand-total/rr2?color=green)](http://www.r-pkg.org/pkg/rr2)

Goal
====

This package provides three R<sup>2</sup>s for statistical models with
correlated errors including classes: ‘lmerMod’ (LMM), ‘glmerMod’ (GLMM), 'gls', ‘phylolm’ (Phylogenetic GLS), and ‘binaryPGLMM/phyloglm/communityPGLMM’
(Phylogenetic Logistic Regression). Detailed technical descriptions can
be found in [Ives 2018](https://doi.org/10.1093/sysbio/syy060).

Installation
============

This package can be installed with:

``` r
install.packages("rr2")

# or install the latest version
# install.packages("devtools")
devtools::install_github("arives/rr2")
```

Package structure
=================

This package has three main functions: `R2.resid()`, `R2.lik()`, and
`R2.pred()`. You can use them individually in the form of, e.g.,
`R2.resid(mod, mod.r)` where `mod` is the full model and `mod.r` is the
reduced model for partial R2s. If you do not include the reduced model
`mod.r`, then the appropriate model with just the intercept is used to
give the total R<sup>2</sup>. When using `R2.resid` and `R2.pred` with
PGLS, you need to include the phylo object containing a phylogenetic
tree, e.g., `R2.resid(mod, mod.r, phy = phy)`.

You can calculate all three R<sup>2</sup>s at the same time with
`R2(mod, mod.r)`. You can also specify which R<sup>2</sup>(s) to
calculate within this function by turning off unwanted methods, e.g.,
`R2(mod, mod.r, resid = FALSE)` or `R2(mod, mod.r, pred = FALSE)`.

This package also has some helper functions such as `inv.logit()`,
`partialR2()`, and `partialR2adj()`.

| Models                           | Available.R2s             |
|:---------------------------------|:--------------------------|
| LM                               | partialR2, partialR2adj   |
| LM                               | R2.pred, R2.resid, R2.lik |
| GLM                              | R2.pred, R2.resid, R2.lik |
| LMM: lmerMod                     | R2.pred, R2.resid, R2.lik |
| GLMM: glmerMod                   | R2.pred, R2.resid, R2.lik |
| PGLS: phylolm                    | R2.pred, R2.resid, R2.lik |
| PGLS: gls                        | R2.pred, R2.resid, R2.lik |
| PGLMM: binaryPGLMM               | R2.pred, R2.resid, ——-    |
| PGLMM: phyloglm                  | ——-, ——–, R2.lik          |
| PGLMM: communityPGLMM (gaussian) | R2.pred, ——–, R2.lik      |
| PGLMM: communityPGLMM (binomial) | R2.pred, ——–, ——-         |

Usage: calculating R<sup>2</sup>s for regression models
=======================================================

First, let’s simulate data that will be used to fit various models.

``` r
# data 
set.seed(123)
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
prob <- rr2::inv.logit(b1 * d$x1 + rep(rnorm(n = p1, sd = sd1), each = nsample)) 
# random intercept u1
d$y_binary <- rbinom(n = n, size = 1, prob = prob)

# PGLS
b1 <- 1.5; signal <- 0.7
phy <- ape::compute.brlen(ape::rtree(n = n), method = "Grafen", power = 1)
phy.x <- ape::compute.brlen(phy, method = "Grafen", power = .0001)
x_trait <- ape::rTraitCont(phy.x, model = "BM", sigma = 1)
e <- signal^0.5 * ape::rTraitCont(phy, model = "BM", sigma = 1) + (1-signal)^0.5 * rnorm(n=n)
d$x_trait <- x_trait[match(names(e), names(x_trait))]
d$y_pgls <- b1 * x_trait + e
rownames(d) <- phy$tip.label    

# Phylogenetic Logistic Regression
b1 <- 1.5; signal <- 2
e <- signal * ape::rTraitCont(phy, model = "BM", sigma = 1)
e <- e[match(phy$tip.label, names(e))]
d$y_phy_binary <- rbinom(n = n, size = 1, prob = rr2::inv.logit(b1 * d$x1 + e))

head(d)
```

    ##              x1          x2 u1 u2 y_re_intercept y_re_slope y_binary
    ## t31 -0.56047565 -0.71040656  1  1       3.053041 -0.2790159        1
    ## t37 -0.23017749  0.25688371  1  2       3.794671  1.7435372        0
    ## t8   1.55870831 -0.24669188  1  3       8.062178 -0.3410566        1
    ## t70  0.07050839 -0.34754260  1  4       3.649759  0.5076822        0
    ## t53  0.12928774 -0.95161857  1  5       2.526704  0.2830316        0
    ## t13  1.71506499 -0.04502772  1  6       7.631604 -8.5551981        0
    ##         x_trait     y_pgls y_phy_binary
    ## t31 -2.07597968 -2.8257102            0
    ## t37 -0.31921893 -0.6918108            0
    ## t8  -0.24097587 -0.3359352            1
    ## t70 -0.08278377 -0.6383157            0
    ## t53 -1.60010819 -1.3718365            0
    ## t13 -1.52297135 -2.0347222            1

Then, let’s fit some models and calculate their R<sup>2</sup>s.

LM
--

``` r
library(rr2)
z.f.lm <- lm(y_re_intercept ~ x1 + x2, data = d)
z.x.lm <- lm(y_re_intercept ~ x1, data = d)
z.0.lm <- lm(y_re_intercept ~ 1, data = d)

R2(mod = z.f.lm, mod.r = z.x.lm)
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.2473776 0.2473776 0.2473776

``` r
partialR2(mod = z.f.lm, mod.r = z.x.lm)
```

    ## [1] 0.2473776

``` r
partialR2adj(mod = z.f.lm, mod.r = z.x.lm)
```

    ## $R2
    ## [1] 0.2473776
    ## 
    ## $R2.adj
    ## [1] 0.4982517

LMM
---

``` r
z.f.lmm <- lme4::lmer(y_re_intercept ~ x1 + x2 + (1 | u1) + (1 | u2), data = d, REML = F)
z.x.lmm <- lme4::lmer(y_re_intercept ~ x1 + (1 | u1) + (1 | u2), data = d, REML = F)
z.v.lmm <- lme4::lmer(y_re_intercept ~ 1 + (1 | u2), data = d, REML = F)
z.0.lmm <- lm(y_re_intercept ~ 1, data = d)

R2(mod = z.f.lmm, mod.r = z.x.lmm)
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.5356524 0.6036311 0.6087728

``` r
R2(mod = z.f.lmm, mod.r = z.v.lmm)
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.7441745 0.8373347 0.8559029

``` r
R2(mod = z.f.lmm, mod.r = z.0.lmm)
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.7762978 0.8767789 0.8991618

``` r
R2(mod = z.f.lmm) # if omit mod.r, default will be the simplest model, such as z.0.lmm here.
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.7762978 0.8767789 0.8991618

GLMM
----

``` r
z.f.glmm <- lme4::glmer(y_binary ~ x1 + (1 | u1), data = d, family = "binomial")
z.x.glmm <- lme4::glmer(y_binary ~ 1 + (1 | u1), data = d, family = "binomial")
z.v.glmm <- glm(y_binary ~ x1, data = d, family = "binomial")

R2(mod = z.f.glmm, mod.r = z.x.glmm)
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.1170588 0.1413694 0.1373521

``` r
R2(mod = z.f.glmm, mod.r = z.v.glmm)
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.1990563 0.3404476 0.3545240

``` r
R2(mod = z.f.glmm)
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.2406380 0.3659939 0.3792381

``` r
# specify sigma2_d for R2.resid()
R2.resid(mod = z.f.glmm, mod.r = z.v.glmm, sigma2_d = "s2w")
```

    ## [1] 0.3404476

``` r
R2.resid(mod = z.f.glmm, mod.r = z.v.glmm, sigma2_d = "NS") 
```

    ## [1] 0.4246935

``` r
R2.resid(mod = z.f.glmm, mod.r = z.v.glmm, sigma2_d = "rNS")
```

    ## [1] 0.4553596

PGLS
----

``` r
z.f.pgls <- phylolm::phylolm(y_pgls ~ x_trait, phy = phy, data = d, model = "lambda")
z.v.lm <- lm(y_pgls ~ x_trait, data = d)

# phy is needed for phylogenetic models' R2.resid and R2.pred
R2(mod = z.f.pgls, mod.r = z.v.lm, phy = phy)
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.2353912 0.3590018 0.3114035

``` r
R2(mod = z.f.pgls, phy = phy)
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.8642865 0.8862266 0.8777782

``` r
# This also works for models fit with nlme::gls()
z.f.gls <- nlme::gls(y_pgls ~ x_trait, data = d, correlation = ape::corPagel(1, phy), method = "ML")
z.x.gls <- nlme::gls(y_pgls ~ 1, data = d, correlation = ape::corPagel(1, phy), method = "ML")
R2(mod = z.f.gls, mod.r = z.v.lm)
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.2353912 0.3590294 0.3114048

``` r
R2(mod = z.f.gls)
```

    ##    R2_lik  R2_resid   R2_pred 
    ## 0.8642865 0.8862315 0.8777784

Phylogenetic Logistic Regression
--------------------------------

*Note*: we modified `ape::binaryPGLMM` to return necessary components
for `rr2::R2()`.

``` r
z.f.plog <- rr2::binaryPGLMM(y_phy_binary ~ x1, data = d, phy = phy)
z.x.plog <- rr2::binaryPGLMM(y_phy_binary ~ 1, data = d, phy = phy)
z.v.plog <- glm(y_phy_binary ~ x1, data = d, family = "binomial")

# R2.lik can't be used with binaryPGLMM because it is not a ML method
R2(mod = z.f.plog, mod.r = z.x.plog)
```

    ## Models of class binaryPGLMM do not have R2.lik method.

    ##  R2_resid   R2_pred 
    ## 0.6115816 0.3344832

``` r
R2(mod = z.f.plog)
```

    ## Models of class binaryPGLMM do not have R2.lik method.

    ##  R2_resid   R2_pred 
    ## 0.8076862 0.5531285

``` r
z.f.plog2 <- phylolm::phyloglm(y_phy_binary ~ x1, data = d, start.alpha = 1, phy = phy)
z.x.plog2 <- phylolm::phyloglm(y_phy_binary ~ 1, data = d, phy = phy, 
                               start.alpha = min(20, z.f.plog2$alpha))
z.v.plog2 <- glm(y_phy_binary ~ x1, data = d, family = "binomial")

# R2.resid and R2.pred do not apply for phyloglm
R2(z.f.plog2, z.x.plog2) 
```

    ## Models of class phyloglm only have R2.lik method.

    ##    R2_lik 
    ## 0.3853273

``` r
# alternate
R2.lik(z.f.plog2, z.x.plog2)
```

    ## [1] 0.3853273

Contributions of predictors
===========================

We can use `rr2::R2()` to calculate partial R<sup>2</sup>s and compare
contributions of different predictors. Here is an example using
`phylolm::phyloglm()`. The same comparisons can be also applied to other
types of models.

``` r
z.f <- phylolm::phyloglm(y_phy_binary ~ x1 + x2, data = d, start.alpha = 1, phy = phy)
z.r1 <- phylolm::phyloglm(y_phy_binary ~ x1, data = d, start.alpha = 1, phy = phy)
z.r2 <- phylolm::phyloglm(y_phy_binary ~ x2, data = d, start.alpha = 1, phy = phy)
# total R2
R2(z.f)
```

    ## Models of class phyloglm only have R2.lik method.

    ##    R2_lik 
    ## 0.4825879

``` r
# contribution of x1
R2(z.f, z.r2)
```

    ## Models of class phyloglm only have R2.lik method.

    ##    R2_lik 
    ## 0.3877523

``` r
# contribution of x2
R2(z.f, z.r1)
```

    ## Models of class phyloglm only have R2.lik method.

    ##      R2_lik 
    ## 0.004148668

It is also possible to estimate the “contribution” of correlation
structrues in the model. For the above example, we can replace the
phylogeny with a star phylogeny and then compare the R<sup>2</sup>s of
the two models.

``` r
# see the first chunk R code for the build of phy.x, a star phylogeny
z.r3 <- phylolm::phyloglm(y_phy_binary ~ x1 + x2, data = d, start.alpha = 1, phy = phy.x)
R2(z.f, z.r3)
```

    ## Models of class phyloglm only have R2.lik method.

    ##    R2_lik 
    ## 0.3124863

Citation
========

Please cite the following papers if you find this package useful:

> -   [Anthony R. Ives. 2018. R2s for Correlated Data: Phylogenetic
>     Models, LMMs, and GLMMs. Systematic Biology,
>     syy060.](https://doi.org/10.1093/sysbio/syy060)  
> -   [Anthony R. Ives and Daijiang Li (2018). rr2: An R package to
>     calculate R^2s for regression models. The Journal of Open Source
>     Software, 3(30), 1028.](https://doi.org/10.21105/joss.01028)

Contributing
============

Contributions are welcome. You can either provide comments and feedback
by filing an issue on Github
[here](https://github.com/arives/rr2/issues) or making pull requests. It
may be easier if you first open an issue outlining what you will do in
the pull request.

Questions about the package can also be posted as issues on Github.
