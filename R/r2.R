#' rr2: An R package to calculate R2s for regression models
#' 
# Copyright (C) 2019 Anthony Ives; Daijiang Li
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' The rr2 package provides methods to calculate R2 for models with correlated errors, 
#' including Phylogenetic GLS, Phylogenetic Logistic Regression, LMMs, GLMM, and PGLMM.
#'  
#' @docType package
#' @name rr2
NULL

#' Calculate R2_lik, R2_resid, and R2_pred
#'
#' This is a wrapper for calculating three R2s -- R2_lik, R2_resid, and R2_pred -- for LMMs and GLMMs, and phylogenetic LMMs (PLMMs) and GLMMs (PGLMMs). Note that the individual functions \code{R2_lik()}, \code{R2_resid()}, and \code{R2_pred()} can be called separately. This is preferrable if you are only interested in one R2; for example, for \code{phylolm()} called from `R2` you need to specify 'phy' (phylo object for the phylogeny), while \code{R2_lik()} does not require this.
#' 
#' Details about the methods are provided under the separate functions for \code{R2_lik()}, \code{R2_resid()}, and \code{R2_pred()}. There are also many worked examples. 
#'   
#' @param mod A regression model with one of the following classes: 'lm', 'glm', lmerMod', glmerMod', 'phylolm', 'gls', 'pglmm', 'pglmm_compare', binaryPGLMM', or 'communityPGLMM'.
#' @param mod.r A reduced model; if not provided, the total R2 will be given by setting 'mod.r' to the model corresponding to 'mod' with the intercept as the only predictor.
#' @param phy The phylogeny for phylogenetic models (as a 'phylo' object), which is not required to be specified for \code{R2_lik()} or non-phylogenetic models.
#' @param sigma2_d Distribution-specific variance \eqn{\sigma^2_d}{sigma2d} (see Details) used in \code{R2_resid()}. For binomial GLMs, GLMMs and PGLMMs with logit link functions, options are c('s2w', 'NS', 'rNS'). For binomial GLMs, GLMMs and PGLMMs with probit link functions, options are c('s2w', 'NS'). Other families use 's2w'.
#' @param lik Whether to calculate R2_lik; default is TRUE.
#' @param resid Whether to calculate R2_resid; default is TRUE.
#' @param pred Whether to calculate R2_pred; default is TRUE.
#' @param ... Additional arguments for `R2_pred()`. `gaussian.pred = "tip_rm"` or `gaussian.pred = "nearest_node"`.
#' @return A vector, with all three R2s by default.
#' @author Daijiang Li and Anthony R. Ives
#' @references Ives A.R. and Li D. 2018. rr2: An R package to calculate R2s for regression models. Journal of Open Source Software. DOI:10.21105/joss.01028
#' 
#' Ives A.R. 2018. R2s for Correlated Data: Phylogenetic Models, LMMs, and GLMMs. Systematic Biology, Volume 68, Issue 2, March 2019, Pages 234-251. DOI:10.1093/sysbio/syy060
#' @seealso MuMIn, lme4, ape, phylolm, phyr, pez
#' @export
#' @examples 
#' library(ape)
#' library(phylolm)
#' library(lme4)
#' library(nlme)
#' library(phyr)
#' 
#' set.seed(12345)
#' p1 <- 10
#' nsample <- 10
#' n <- p1 * nsample
#' 
#' d <- data.frame(x1 = 0, x2 = 0, u1 = rep(1:p1, each = nsample),
#'                 u2 = rep(1:p1, times = nsample))
#' d$u1 <- as.factor(d$u1)
#' d$u2 <- as.factor(d$u2)
#' 
#' b1 <- 1
#' b2 <- -1
#' sd1 <- 1.5
#' 
#' d$x1 <- rnorm(n = n)
#' d$x2 <- rnorm(n = n)
#' d$y.lmm <- b1 * d$x1 + b2 * d$x2 + 
#'   rep(rnorm(n = p1, sd = sd1), each = nsample) +
#'   rep(rnorm(n = p1, sd = sd1), times = nsample) + 
#'   rnorm(n = n)
#' 
#' prob <- inv.logit(b1 * d$x1 + rep(rnorm(n = p1, sd = sd1), each = nsample))
#' d$y.glmm <- rbinom(n = n, size = 1, prob = prob)
#' 
#' # LMM with two fixed and two random effects ----
#' z.f <- lmer(y.lmm ~ x1 + x2 + (1 | u1) + (1 | u2), data = d, REML = FALSE)
#' z.x <- lmer(y.lmm ~ x1 + (1 | u1) + (1 | u2), data = d, REML = FALSE)
#' z.v <- lmer(y.lmm ~ 1 + (1 | u2), data = d, REML = FALSE)
#' z.0 <- lm(y.lmm ~ 1, data = d)
#' 
#' R2(z.f, z.x)
#' R2(z.f, z.v)
#' R2(z.f)
#' 
#' # GLMM with one fixed and one random effect ----
#' z.f <- glmer(y.glmm ~ x1 + (1 | u1), data = d, family = 'binomial')
#' z.x <- glmer(y.glmm ~ 1 + (1 | u1), data = d, family = 'binomial')
#' z.v <- glm(y.glmm ~ x1, data = d, family = 'binomial')
#' 
#' R2(z.f, z.x)
#' R2(z.f, z.v)
#' R2(z.f)
#' 
#' # These give different results for R2_resid.
#' R2(z.f, sigma2_d = 's2w')
#' R2(z.f, sigma2_d = 'NS')
#' R2(z.f, sigma2_d = 'rNS')
#' 
#' # GLS {nlme} with one fixed effect and autocorrelated errors among 6 groups ----
#' nT <- 10
#' nseries <- 6
#' n <- nT * nseries
#' 
#' d <- data.frame(x = 0, y = 0, u = rep(1:nseries, each = nT), e = rnorm(1))
#' d$u <- as.factor(d$u)
#' d$x <- rnorm(n = n)
#' ar1 <- .5
#' for(t in 2:n) d$e[t] <- ar1*d$e[t-1] + rnorm(1)
#' 
#' b1 <- 1
#' d$y <- b1 * d$x + d$e
#' 
#' z.f <- gls(y ~ x + u, correlation = corAR1(form = ~1 | u), data = d)
#' z.x <- gls(y ~ 1, correlation = corAR1(form = ~1 | u), data = d)
#' z.ar <- lm(y ~ x + u, data = d)
#' 
#' R2(z.f, z.x)
#' R2(z.f, z.ar)
#' R2(z.f)
#' 
#' # PGLS with a single fixed effect ----
#' n <- 100
#' d <- data.frame(x = rep(0, n), y = 0)
#' 
#' b1 <- 1.5
#' signal <- 0.7
#' 
#' phy <- compute.brlen(rtree(n = n), method = 'Grafen', power = 1)
#' phy.x <- compute.brlen(phy, method = 'Grafen', power = .0001)
#' 
#' # Generate random data
#' x <- rTraitCont(phy.x, model = 'BM', sigma = 1)
#' e <- signal ^ 0.5 * rTraitCont(phy, model = 'BM', sigma = 1) +
#'   (1 - signal) ^ 0.5 * rnorm(n = n)
#' d$x <- x[match(names(e), names(x))]
#' d$y <- b1 * x + e
#' rownames(d) <- phy$tip.label
#' d$sp <- phy$tip.label
#' 
#' # Fit with phylolm() in {phylolm}
#' z.f <- phylolm(y ~ x, phy = phy, data = d, model = 'lambda')
#' z.x <- phylolm(y ~ 1, phy = phy, data = d, model = 'lambda')
#' z.v <- lm(y ~ x, data = d)
#' 
#' R2(z.f, z.x, phy = phy)
#' R2(z.f, z.v, phy = phy)
#' R2(z.f, phy = phy)
#' 
#' # These data can also be fit with pglmm_compare in {phyr}
#' # Note that pglmm_compare will be renamed to pglmm_compare in the next version
#' z.f <- pglmm_compare(y ~ x, data = d, phy = phy, REML=FALSE)
#' z.x <- pglmm_compare(y ~ 1, data = d, phy = phy, REML=FALSE)
#' z.v <- glm(y ~ x, data = d)
#' 
#' R2(z.f, z.x)
#' R2(z.f, z.v)
#' R2(z.f)
#' 
#' # This also works for models fit with gls() in {nlme}
#' z.f <- gls(y ~ x, data = d, correlation = corPagel(1, phy, form = ~sp), method = "ML")
#' z.x <- gls(y ~ 1, data = d, correlation = corPagel(1, phy, form = ~sp), method = "ML")
#' z.v <- lm(y ~ x, data = d)
#' 
#' R2(z.f, z.x)
#' R2(z.f, z.v)
#' R2(z.f)
#' 
#' # But note that you need to define weights for gls() with non-ultrametric trees;
#' # if not, you will get a error from R2_resid,  "Matrix is not block-diagonal"
#' 
#' phy.nu <- rtree(n = n)
#' # Generate random data
#' e <- signal ^ 0.5 * rTraitCont(phy.nu, model = 'BM', sigma = 1) +
#'   (1 - signal) ^ 0.5 * rnorm(n = n)
#' d$x <- x[match(names(e), names(x))]
#' d$y <- b1 * x + e
#' rownames(d) <- phy.nu$tip.label
#' d$sp <- phy.nu$tip.label
#' 
#' weights <- diag(vcv.phylo(phy.nu))
#' z.f <- gls(y ~ x, data = d,
#'            correlation = corPagel(1, phy.nu, form = ~sp),
#'            weights = varFixed(~weights), method = "ML")
#' z.x <- gls(y ~ 1, data = d,
#'            correlation = corPagel(1, phy.nu, form = ~sp),
#'            weights = varFixed(~weights), method = "ML")
#' z.v <- lm(y ~ x, weights = 1/weights, data = d)
#' 
#' R2(z.f, z.x)
#' R2(z.f, z.v)
#' R2(z.f)
#' 
#' # PGLMM with one fixed effect ----
#' n <- 100
#' b1 <- 1.5
#' signal <- 2
#' 
#' phy <- compute.brlen(rtree(n = n), method = 'Grafen', power = 1)
#' phy.x <- compute.brlen(phy, method = 'Grafen', power = .0001)
#' 
#' # Generate random data
#' x <- rnorm(n)
#' d <- data.frame(x = x, y = 0)
#' 
#' e <- signal * rTraitCont(phy, model = 'BM', sigma = 1)
#' e <- e[match(phy$tip.label, names(e))]
#' 
#' d$y <- rbinom(n = n, size = 1, prob = inv.logit(b1 * d$x + e))
#' rownames(d) <- phy$tip.label
#' # Use the function phyloglm() from the phylolm package.
#' z.f <- phyloglm(y ~ x, data = d, start.alpha = 1, phy = phy)
#' z.x <- phyloglm(y ~ 1, data = d, phy = phy, start.alpha = min(20, z.f$alpha))
#' z.v <- glm(y ~ x, data = d, family = 'binomial')
#' 
#' R2(z.f, z.x)
#' R2(z.f, z.v)
#' R2(z.f)
#' 
#' # Use the function pglmm_compare() from the phyr package. Note that this is a 
#' # different model from phyloglm()
#' z.f <- pglmm_compare(y ~ x, data = d, family = 'binomial', phy = phy, REML = FALSE)
#' z.x <- pglmm_compare(y ~ 1, data = d, family = 'binomial', phy = phy, REML = FALSE)
#' z.v <- glm(y ~ x, data = d, family = 'binomial')
#' 
#' R2(z.f, z.x)
#' R2(z.f, z.v)
#' R2(z.f)
#' 
#' # A community example of pglmm {phyr} ----
#' library(mvtnorm)
#' nspp <- 6
#' nsite <- 4
#' 
#' # Simulate a phylogeny that has a lot of phylogenetic signal (power = 1.3)
#' phy <- compute.brlen(rtree(n = nspp), method = "Grafen", power = 1.3)
#' 
#' # Simulate species means
#' sd.sp <- 1
#' mean.sp <- rTraitCont(phy, model = "BM", sigma=sd.sp ^ 2)
#' 
#' # Replicate values of mean.sp over sites
#' Y.sp <- rep(mean.sp, times = nsite)
#' 
#' # Simulate site means
#' sd.site <- 1
#' mean.site <- rnorm(nsite, sd = sd.site)
#' 
#' # Replicate values of mean.site over sp
#' Y.site <- rep(mean.site, each = nspp)
#' 
#' # Compute a covariance matrix for phylogenetic attraction
#' sd.attract <- 1
#' Vphy <- vcv(phy)
#' 
#' # Standardize the phylogenetic covariance matrix to have determinant = 1.
#' # (For an explanation of this standardization, see subsection 4.3.1 in Ives (2018))
#' Vphy <- Vphy/(det(Vphy)^(1/nspp))
#' 
#' # Construct the overall covariance matrix for phylogenetic attraction.
#' # (For an explanation of Kronecker products, see subsection 4.3.1 in the book)
#' V <- kronecker(diag(nrow = nsite, ncol = nsite), Vphy)
#' Y.attract <- array(t(rmvnorm(n = 1, sigma = sd.attract^2*V)))
#' 
#' # Simulate residual errors
#' sd.e <- 1
#' Y.e <- rnorm(nspp * nsite, sd = sd.e)
#' 
#' # Construct the dataset
#' d <- data.frame(sp = rep(phy$tip.label, times = nsite), site = rep(1:nsite, each = nspp))
#' 
#' # Simulate abundance data
#' d$Y <- Y.sp + Y.site + Y.attract + Y.e
#' 
#' # Full and reduced models
#' z.f <- pglmm(Y ~ 1 + (1|sp__) + (1|site) + (1|sp__@site),
#'              data = d, cov_ranef = list(sp = phy), REML = FALSE)
#' z.nested <- pglmm(Y ~ 1 + (1|sp__) + (1|site),
#'                   data = d, cov_ranef = list(sp = phy), REML = FALSE)
#' z.sp <- pglmm(Y ~ 1 + (1|sp) + (1|site),
#'               data = d, cov_ranef = list(sp = phy), REML = FALSE)
#' 
#' R2(z.f, z.nested)
#' R2(z.nested, z.sp)
#' R2(z.f)

R2 <- function(mod = NULL, mod.r = NULL, phy = NULL, sigma2_d = c("s2w", "NS", "rNS"), 
    lik = TRUE, resid = TRUE, pred = TRUE, ...) {
    
    if (all(!lik, !resid, !pred)) 
        stop("Specify at least one of 'lik', 'resid', or 'pred' for R2")
    
    sigma2_d <- match.arg(sigma2_d)
    
    # phyloglm only has a R2_lik method.
    if (any(class(mod) %in% "phyloglm")) {
        resid <- FALSE
        pred <- FALSE
        message("Models of class phyloglm only have a R2_lik method.")
    }
    
    # binaryPGLMM does not have a R2_lik method.
    if (any(class(mod) %in% "binaryPGLMM") & lik == TRUE) {
        lik <- FALSE
        message("Models of class binaryPGLMM do not have a R2_lik method.")
    }
    
    # pglmm does not have a R2_resid method
    if (any(class(mod) %in% c("communityPGLMM", "pglmm")) & (resid == TRUE | lik == TRUE)) {
      message("Models of class pglmm do not have a R2_resid method.")
      resid <- FALSE
      if (mod$bayes == TRUE) {
        message("Models of class pglmm with bayes = TRUE do not have a R2_lik method.")
        lik <- FALSE
      }
    }
    
    # gls does not have a R2_resid method
    if (any(class(mod) %in% c("gls")) && 
        !any(class(mod$modelStruct$corStruct) %in% 
             c("corBrownian", "corMartins", "corGrafen", "corPagel", "corBlomberg")) && 
        (resid == TRUE)) {
      message("Models of class gls that are not phylogenetic do not have a R2_resid method.")
      resid <- FALSE
    }
    
    out <- rep(NA, 3)
    names(out) <- c("R2_lik", "R2_resid", "R2_pred")
    if (lik) 
        out[1] <- R2_lik(mod, mod.r)
    if (resid) 
        out[2] <- R2_resid(mod, mod.r, sigma2_d, phy = phy)
    if (pred) 
        out[3] <- R2_pred(mod, mod.r, phy = phy, ...)
    
    out <- out[!is.na(out)]  # remove R2s not calculated
    
    return(out)
}

