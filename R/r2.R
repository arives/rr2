#' Calculate R2.lik, R2.resid, and R2.pred
#'
#' This is a wrapper for calculating all three R2s -- R2.lik, R2.resid, and R2.pred -- for LMM, GLMM, PGLM, and PGLMM. Note that the individual functions R2.lik(), R2.resid(), and R2.pred() can be called separately. This is preferrable if you are only  interested in one R2; for example, for phylolm() called from `R2` you need to specify 'phy' (phylo object for the phylogeny), while R2.lik() does not require this.
#' 
#' Details about the methods are provided under the separate functions for `R2.lik`, `R2.resid`, and `R2.pred`. There are also many worked examples. 
#'   
#' @param mod A regression model with the following class: 'lm', 'glm', lmerMod', glmerMod', 'phylolm', 'binaryPGLMM', or 'communityPGLMM'.
#' @param mod.r A reduced model; if not provided, the total R2 will be given by setting 'mod.r' to the model corresponding to 'mod' with intercept as the only predictor.
#' @param phy The phylogeny for phylogenetic models (as a 'phylo' object), which is not required to be specified for R2.lik().
#' @param sigma2_d Distribution-specific variance σ2d (see Details) used in R2.resid(). For binomial GLMs, GLMMs and PGLMMs with logit link functions, options are `c("s2w", "NS", "rNS")`. For binomial GLMs, GLMMs and PGLMMs with probit link functions, options are `c("s2w", "NS")`. Other families use `"s2w"`.
#' @param lik Whether to calculate R2.lik; default is TRUE.
#' @param resid Whether to calculate R2.resid; default is TRUE.
#' @param pred Whether to calculate R2.pred; default is TRUE.
#' @return An array, with all three R2s by default.
#' @author Daijiang Li and Anthony R. Ives
#' @references Ives A. in press. R2s for Correlated Data: Phylogenetic Models, LMMs, and GLMMs. Systematic Biology.
#' @export
#' @examples library(ape)
#' library(phylolm)
#' library(lme4)
#' 
#' #################
#' # LMM with two fixed and two random effects 
#' p1 <- 10
#' nsample <- 10
#' n <- p1 * nsample
#' 
#' d <- data.frame(x1=0, x2=0, y=0, u1=rep(1:p1, each=nsample), u2=rep(1:p1, times=nsample))
#' d$u1 <- as.factor(d$u1)
#' d$u2 <- as.factor(d$u2)
#' 
#' b1 <- 1
#' b2 <- -1
#' sd1 <- 1.5
#' 
#' d$x1 <- rnorm(n=n)
#' d$x2 <- rnorm(n=n)
#' d$y <- b1 * d$x1 + b2 * d$x2 + rep(rnorm(n=p1, sd=sd1), each=nsample) + 
#'        rep(rnorm(n=p1, sd=sd1), times=nsample) + rnorm(n=n)
#' 
#' z.f <- lmer(y ~ x1 + x2 + (1 | u1) + (1 | u2), data=d, REML = FALSE)
#' z.x <- lmer(y ~ x1 + (1 | u1) + (1 | u2), data=d, REML = FALSE)
#' z.v <- lmer(y ~ 1 + (1 | u2), data=d, REML = FALSE)
#' z.0 <- lm(y ~ 1, data=d)
#' 
#' R2(z.f, z.x)
#' R2(z.f, z.v)
#' R2(z.f)
#' 
#' #################
#' # GLMM with one fixed and one random effect
#'
#' p1 <- 10
#' nsample <- 10
#' n <- p1 * nsample
#' 
#' d <- data.frame(x=0, y=0, u=rep(1:p1, each=nsample))
#' d$u <- as.factor(d$u)
#' 
#' b1 <- 1
#' sd1 <- 1.5
#' 
#' d$x <- rnorm(n=n)
#' prob <- inv.logit(b1 * d$x + rep(rnorm(n=p1, sd=sd1), each=nsample))
#' d$y <- rbinom(n=n, size=1, prob=prob)
#' 
#' z.f <- glmer(y ~ x + (1 | u), data=d, family="binomial")
#' z.x <- glmer(y ~ 1 + (1 | u), data=d, family="binomial")
#' z.v <- glm(y ~ x, data=d, family="binomial")
#' 
#' R2(z.f, z.x)
#' R2(z.f, z.v)
#' R2(z.f)
#' 
#' # These give different results for R2.resid
#' R2(z.f, sigma2_d="s2w")
#' R2(z.f, sigma2_d="NS")
#' R2(z.f, sigma2_d="rNS")
#' 
#' #################
#' # PGLS with a single fixed effect
#' 
#' n <- 100
#' d <- data.frame(x=array(0, dim=n), y=0)
#' 
#' b1 <- 1.5
#' signal <- 0.7
#' 
#' phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
#' phy.x <- compute.brlen(phy, method = "Grafen", power = .0001)
#' 
#' # Generate random data
#' x <- rTraitCont(phy.x, model = "BM", sigma = 1)
#' e <- signal^0.5 * rTraitCont(phy, model = "BM", sigma = 1) + (1-signal)^0.5 * rnorm(n=n)
#' d$x <- x[match(names(e), names(x))]
#' d$y <- b1 * x + e
#' rownames(d) <- phy$tip.label	
#' 
#' z.x <- phylolm(y ~ 1, phy=phy, data=d, model="lambda")
#' lam.x <- round(z.x$optpar, digits=4)
#' z.f <- phylolm(y ~ x, phy=phy, data=d, model="lambda")
#' z.v <- lm(y ~ x, data=d)
#' 
#' R2(z.f, z.x, phy = phy)
#' R2(z.f, z.v, phy = phy)
#' R2(z.f, phy = phy)
#' 
#' #################
#' # PGLMM with one fixed effect
#' 
#' n <- 100
#' b1 <- 1.5
#' signal <- 2
#' 
#' phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
#' phy.x <- compute.brlen(phy, method = "Grafen", power = .0001)
#' 
#' # Generate random data
#' x <- rnorm(n)
#' d <- data.frame(x=x, y=0)
#' 
#' e <- signal * rTraitCont(phy, model = "BM", sigma = 1)
#' e <- e[match(phy$tip.label, names(e))]
#' 
#' d$y <- rbinom(n=n, size=1, prob=inv.logit(b1 * d$x + e))
#' rownames(d) <- phy$tip.label	
#' 
#' # Use the function binaryPGLMM() from the rr2 package rather than ape.
#' z.f <- rr2::binaryPGLMM(y ~ x, data=d, phy=phy)
#' z.x <- rr2::binaryPGLMM(y ~ 1, data=d, phy=phy)
#' z.v <- glm(y ~ x, data=d, family="binomial")
#' 
#' # R2.lik is not produced, because binaryPGLMM() does not generate a likelihood.
#' R2(z.f, z.x, phy = phy)
#' R2(z.f, z.v, phy = phy)
#' R2(z.f, phy = phy)
#'
R2 <- function(mod = NULL, mod.r = NULL, phy = NULL, sigma2_d = c("s2w", "NS", "rNS"), 
               lik = TRUE, resid = TRUE, pred = TRUE) {
  
  if(all(!lik, !resid, !pred)) 
    stop("Specify at least one of 'lik', 'resid', or 'pred' for R2")
  
  sigma2_d <- match.arg(sigma2_d)
  
  # phyloglm only have R2.lik method
  if (any(class(mod) %in% "phyloglm")) {
    resid <- FALSE
    pred <- FALSE
    message("Models of class phyloglm only have R2.lik method")
  }
  
  # gaussian communityPGLMM only have R2.lik method
  if (any(class(mod) %in% "communityPGLMM")) {
    if(mod$family == "gaussian"){
      resid <- FALSE
      message("Models of class communityPGLMM (gaussian) do not have R2.resid method")
    }
  }
  
  # binaryPGLMM does not have R2.lik method
  if (any(class(mod) %in% "binaryPGLMM") & lik == TRUE) {
    lik <- FALSE
    message("Models of class binaryPGLMM do not have R2.lik method")
  }
  
  # binary communityPGLMM only have R2.pred method at this moment
  if (any(class(mod) %in% "communityPGLMM") & (lik == TRUE | resid == TRUE)) {
    if(mod$family == "binomial"){
      resid <- FALSE
      lik <- FALSE
      message("Models of class communityPGLMM (binomial) only have R2.pred method")
    }
  }
  
  # phylolm requires phy object except for R2.lik
  if (any(class(mod) %in% "phylolm") & is.null(phy)) 
    stop("Phy object is required for models with class phylolm")
  
  out <- rep(NA, 3)
  names(out) <- c("R2_lik", "R2_resid", "R2_pred")
  if (lik) 
    out[1] <- R2.lik(mod, mod.r)
  if (resid) 
    out[2] <- R2.resid(mod, mod.r, phy, sigma2_d)
  if (pred) 
    out[3] <- R2.pred(mod, mod.r, phy)
  
  out <- out[!is.na(out)] # remove R2s not calculated
  
  return(out)
}
