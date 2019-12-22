#' Calculate R2_lik
#'
#' Calculate partial and total R2s for LMM, GLMM, PGLS, and PGLMM using R2_lik, an R2 based on the likelihood of observing the data.
#' 
#' @param mod A regression model with one of the following classes: 'lm', 'glm', 'lmerMod', 'glmerMod', 'phylolm', 'phyloglm', 'gls', 'pglmm', pglmm_compare' or 'communityPGLMM'.
#' @param mod.r A reduced model; if not provided, the total R2 will be given by setting 'mod.r' to the model corresponding to 'mod' with the intercept as the only predictor.
#' @return R2_lik value.
#' @export
#'
#' @details  \code{R2_lik()} works with classes 'lm', 'glm', 'lmerMod', 'glmerMod', 'phylolm', 'phyloglm', and 'communityPGLMM' (family =  'gaussian' only). It is implemented as
#' 
#' \deqn{partial R2 = 1 - exp(-2/n * (logLik(mod.f) - logLik(mod.r)))}
#' 
#' where 'mod.f' and 'mod.r' are the full and reduced models, respectively. The total R2 is given when 'mod.r' is the model corresponding to mod.f that contains only the intercept. For GLMMs and PGLMMs, \code{R2_lik()} is standardized to have a maximum of one following Nagelkerke (1991).
#' 
#' Note that \code{phyloglm()} can have difficulties in finding solutions when there is no phylogenetic signal. Therefore, when alphaWarn == 2, indicating no phylogenetic signal, the model is refit with the corresponding GLM.
#' 
#' \code{R2_lik()} is also computed for LMMs and GLMMs in the {MuMIn} package.
#' 
#' @seealso MuMIn, lme4, ape, phylolm, phyr, pez
#' @author Anthony R. Ives
#' @references Ives A.R. and Li D. 2018. rr2: An R package to calculate R2s for regression models. Journal of Open Source Software. DOI:10.21105/joss.01028
#' 
#' Ives A.R. 2018. R2s for Correlated Data: Phylogenetic Models, LMMs, and GLMMs. Systematic Biology, Volume 68, Issue 2, March 2019, Pages 234-251. DOI:10.1093/sysbio/syy060
#' 
#' Nagelkerke 1991. A note on a general definition of the coefficient of determination. Biometrika 78:691–692.
#' 
#' @examples 
#' library(ape)
#' library(phylolm)
#' library(lme4)
#' library(nlme)
#' library(phyr)
#' 
#' #################
#' # LMM with two fixed and two random effects 
#' p1 <- 10
#' nsample <- 10
#' n <- p1 * nsample
#' 
#' d <- data.frame(x1 = 0, x2 = 0, y = 0, u1 = rep(1:p1, each = nsample), 
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
#' d$y <- b1 * d$x1 + b2 * d$x2 + rep(rnorm(n = p1, sd = sd1), each = nsample) + 
#'        rep(rnorm(n = p1, sd = sd1), times = nsample) + rnorm(n = n)
#' 
#' z.f <- lmer(y ~ x1 + x2 + (1 | u1) + (1 | u2), data = d, REML = FALSE)
#' z.x <- lmer(y ~ x1 + (1 | u1) + (1 | u2), data = d, REML = FALSE)
#' z.v <- lmer(y ~ 1 + (1 | u2), data = d, REML = FALSE)
#' z.0 <- lm(y ~ 1, data = d)
#' 
#' R2_lik(z.f, z.x)
#' R2_lik(z.f, z.v)
#' R2_lik(z.f)
#' 
#' # These give the same results.
#' R2_lik(z.f, z.0)
#' R2_lik(z.f)
#' 
#' #################
#' # GLMM with one fixed and one random effect
#'
#' p1 <- 10
#' nsample <- 10
#' n <- p1 * nsample
#' 
#' d <- data.frame(x = 0, y = 0, u = rep(1:p1, each = nsample))
#' d$u <- as.factor(d$u)
#' 
#' b1 <- 1
#' sd1 <- 1.5
#' 
#' d$x <- rnorm(n = n)
#' prob <- inv.logit(b1 * d$x + rep(rnorm(n = p1, sd = sd1), each = nsample))
#' d$y <- rbinom(n = n, size = 1, prob = prob)
#' 
#' z.f <- glmer(y ~ x + (1 | u), data = d, family = 'binomial')
#' z.x <- glmer(y ~ 1 + (1 | u), data = d, family = 'binomial')
#' z.v <- glm(y ~ x, data = d, family = 'binomial')
#' 
#' R2_lik(z.f, z.x)
#' R2_lik(z.f, z.v)
#' R2_lik(z.f)
#' 
#' #################
#' # PGLS with a single fixed effect
#' 
#' n <- 100
#' d <- data.frame(x = array(0, dim = n), y = 0)
#' 
#' b1 <- 1.5
#' signal <- 0.7
#' 
#' phy <- compute.brlen(rtree(n = n), method = 'Grafen', power = 1)
#' phy.x <- compute.brlen(phy, method = 'Grafen', power = .0001)
#' 
#' # Generate random data
#' x <- rTraitCont(phy.x, model = 'BM', sigma = 1)
#' e <- signal^0.5 * rTraitCont(phy, model = 'BM', sigma = 1) + (1-signal)^0.5 * rnorm(n = n)
#' d$x <- x[match(names(e), names(x))]
#' d$y <- b1 * x + e
#' rownames(d) <- phy$tip.label
#' 
#' z.x <- phylolm(y ~ 1, phy = phy, data = d, model = 'lambda')
#' lam.x <- round(z.x$optpar, digits = 4)
#' z.f <- phylolm(y ~ x, phy = phy, data = d, model = 'lambda')
#' z.v <- lm(y ~ x, data = d)
#' 
#' R2_lik(z.f, z.x)
#' R2_lik(z.f, z.v)
#' R2_lik(z.f)
#' 
#' # This also works for models fit with gls() in {nlme}
#' z.x <- gls(y ~ 1, data = d, correlation = corPagel(1, phy), method = "ML")
#' z.f <- gls(y ~ x, data = d, correlation = corPagel(1, phy), method = "ML")
#' z.v <- lm(y ~ x, data = d)
#' R2_lik(z.f, z.x)
#' R2_lik(z.f, z.v)
#' R2_lik(z.f)
#' 
#' #################
#' # PGLMM with one fixed effect
#' 
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
#' 
#' z.f <- phyloglm(y ~ x, data = d, start.alpha = 1, phy = phy)
#' z.x <- phyloglm(y ~ 1, data = d, phy = phy, start.alpha = min(20,z.f$alpha))
#' z.v <- glm(y ~ x, data = d, family = 'binomial')
#' 
#' R2_lik(z.f, z.x)
#' R2_lik(z.f, z.v)
#' R2_lik(z.f)
#' 
#' # These data can also be fit with pglmm_compare(), although note that this is a different model from phyloglm()
#' z.f <- pglmm_compare(y ~ x, data = d, family = "binomial", phy = phy, REML=F)
#' z.x <- pglmm_compare(y ~ 1, data = d, family = "binomial", phy = phy, REML=F)
#' z.v <- glm(y ~ x, data = d, family = "binomial")
#' 
#' R2_lik(z.f, z.x)
#' R2_lik(z.f, z.v)
#' R2_lik(z.f)

R2_lik <- function(mod = NULL, mod.r = NULL) {
    if (class(mod)[1] == "merModLmerTest") 
        class(mod) <- "lmerMod"
    
    if (!is.element(class(mod)[1], c("lm", "glm", "lmerMod", "glmerMod", "phylolm", 
        "phyloglm", "gls", "pglmm", "pglmm_compare", "communityPGLMM"))) {
        stop("mod must be class one of classes lm, glm, lmerMod, glmerMod, phylolm, phyloglm, gls, pglmm, pglmm_compare, communityPGLMM.")
    }
    
    if (class(mod)[1] == "lm") {
        if (!is.object(mod.r)) {
            y <- model.frame(mod)[, 1]
            mod.r <- lm(y ~ 1)
        }
        if (!is.element(class(mod.r)[1], c("lm"))) {
            stop("mod.r must be class lm.")
        }
        return(R2_lik.lm(mod, mod.r)[1])
    }
    
    if (class(mod)[1] == "glm") {
        if (!is.object(mod.r)) {
            y <- model.frame(mod)[, 1]
            mod.r <- glm(y ~ 1, family = family(mod)[[1]])
        }
        if (!is.element(class(mod.r)[1], c("glm"))) {
            stop("mod.r must be class glm.")
        }
        if (family(mod)[[1]] != family(mod.r)[[1]]) {
            stop("Sorry, but mod and mod.r must be from the same family of distributions.")
        }
        return(R2_lik.glm(mod, mod.r)[1])
    }
    
    if (class(mod)[1] == "lmerMod") {
        if (!is.object(mod.r)) {
            y <- model.frame(mod)[, 1]
            mod.r <- lm(y ~ 1)
        }
        if (class(mod.r)[1] == "merModLmerTest") 
            class(mod.r) <- "lmerMod"
        if (!is.element(class(mod.r)[1], c("lmerMod", "lm"))) {
            stop("mod.r must be class lmerMod or lm.")
        }
        if (lme4::isREML(mod)) {
            mod <- update(mod, REML = F)
            warning("mod updated with REML = F")
        }
        if (class(mod.r)[1] == "lmerMod" && lme4::isREML(mod.r)) {
            mod.r <- update(mod.r, REML = F)
            warning("mod.r updated with REML = F")
        }
        return(R2_lik.lmerMod(mod, mod.r)[1])
    }
    
    if (class(mod)[1] == "glmerMod") {
        if (!is.object(mod.r)) {
            y <- model.frame(mod)[, 1]
            mod.r <- glm(y ~ 1, family = family(mod)[[1]])
        }
        if (!is.element(class(mod.r)[1], c("glmerMod", "glm"))) {
            stop("mod.r must be class glmerMod or glm.")
        }
        if (family(mod)[[1]] != family(mod.r)[[1]]) {
            stop("Sorry, but mod and mod.r must be from the same family of distributions.")
        }
        return(R2_lik.glmerMod(mod, mod.r)[1])
    }
    
    if (class(mod)[1] == "phylolm") {
        if (!is.object(mod.r)) {
            y <- mod$y
            mod.r <- lm(y ~ 1)
        }
        if (!is.element(class(mod.r)[1], c("phylolm", "lm"))) {
            stop("mod.r must be class phylolm or lm.")
        }
        return(R2_lik.phylolm(mod, mod.r))
    }
    
    if (class(mod)[1] == "phyloglm") {
        if (!is.object(mod.r)) {
            y <- mod$y
            mod.r <- glm(y ~ 1, family = "binomial")
        }
        if (!is.element(class(mod.r)[1], c("phyloglm", "glm"))) {
            stop("mod.r must be class phyloglm or glm.")
        }
        return(R2_lik.phyloglm(mod, mod.r)[1])
    }
    
    if (class(mod)[1] == "gls") {
      if (!is.object(mod.r)) {
        y <- as.numeric(fitted(mod) + resid(mod))
        mod.r <- lm(y ~ 1)
      }
      if (!is.element(class(mod.r)[1], c("gls", "lm"))) {
        stop("mod.r must be class gls or lm.")
      }
      return(R2_lik.gls(mod, mod.r))
    }

    if (class(mod)[1] %in% c("pglmm", "pglmm_compare")) {
      if (mod$REML == TRUE) 
        stop("mod was fitted with REML, please set it to FALSE and re-fit it")
      
      if (!is.object(mod.r)) {
        y <- mod$Y
        mod.r <- glm(y ~ 1, family = mod$family)
      }
      if (!is.element(class(mod.r)[1], c("pglmm", "pglmm_compare", "glm"))) {
        stop("mod.r must be class pglmm, pglmm_compare or glm.")
      }
      return(R2_lik.pglmm(mod, mod.r))
    }

    if (class(mod)[1] == "communityPGLMM") {
       if (mod$family == "binomial") 
         stop("Binary communityPGLMMs do not have log likelihood,
                  If you are interested in LRT of random terms, use
                  phyr::communityPGLMM.binary.LRT()")
      if (mod$REML == TRUE) 
         stop("mod was fitted with REML, please set it to FALSE and re-fit it")
        
      if (!is.object(mod.r)) {
            y <- mod$Y
            mod.r <- lm(y ~ 1)
      }
      if (!is.element(class(mod.r)[1], c("communityPGLMM", "lm"))) {
          stop("mod.r must be class communityPGLMM or lm.")
      }
      return(R2_lik.communityPGLMM(mod, mod.r))
    }
}

R2_lik.lm <- function(mod = NULL, mod.r = NULL) {
    X <- model.matrix(mod)
    n <- dim(X)[1]
    
    R2_lik <- 1 - exp(-2/n * (logLik(mod)[[1]] - logLik(mod.r)[[1]]))
    return(R2_lik)
}

R2_lik.glm <- R2_lik.lm

R2_lik.lmerMod <- R2_lik.lm

R2_lik.glmerMod <- function(mod = NULL, mod.r = NULL) {
    X <- model.matrix(mod)
    n <- dim(X)[1]
    R2_lik <- (1 - exp(-2/n * (logLik(mod) - logLik(mod.r))))/(1 - exp(2/n * logLik(mod.r)))
    return(R2_lik)
}

R2_lik.phylolm <- function(mod = NULL, mod.r = NULL) {
    X <- mod$X
    n <- dim(X)[1]
    R2_lik <- 1 - exp(-2/n * (logLik(mod)[[1]] - logLik(mod.r)[[1]]))
    return(R2_lik)
}

R2_lik.phyloglm <- function(mod = NULL, mod.r = NULL) {
    
    y <- mod$y
    X <- mod$X
    n <- dim(X)[1]
    
    if (mod$alphaWarn == 2) {
      LL <- logLik(glm(y ~ 0 + X, family = "binomial"))
      warning("In mod, alphaWarn = 2, so model refit with glm()")
    } else {
      LL <- mod$logLik
    }
    if (class(mod.r)[1] == "phyloglm") {
        if (mod.r$alphaWarn == 2) {
          X.r <- mod.r$X
          LL.r <- logLik(glm(y ~ 0 + X.r, family = "binomial"))
          warning("In mod.r, alphaWarn = 2, so model refit with glm()")
        } else {
          LL.r <- mod.r$logLik
        }
    } else {
        LL.r <- logLik(mod.r)
    }
    
    R2_lik <- (1 - exp(-2/n * (LL - LL.r)))/(1 - exp(2/n * LL.r))
    
    return(R2_lik)
}

R2_lik.gls <- function(mod = NULL, mod.r = NULL) {
  n <- mod$dims$N
  R2_lik <- 1 - exp(-2/n * (logLik(mod)[[1]] - logLik(mod.r)[[1]]))
  return(R2_lik)
}

R2_lik.pglmm <- function(mod = NULL, mod.r = NULL) {
  n <- nrow(mod$X)
  if (any(class(mod.r) %in% c("lm", "glm"))) {
    ll.r <- logLik(mod.r)[[1]]
  } else {
    ll.r <- mod.r$logLik
  }
  R2_lik <- 1 - exp(-2/n * (mod$logLik - ll.r))
  return(R2_lik)
}

R2_lik.communityPGLMM <- function(mod = NULL, mod.r = NULL) {
    n <- nrow(mod$X)
    if (class(mod.r) == "lm") {
        ll.r <- logLik(mod.r)[[1]]
    } else {
        ll.r <- mod.r$logLik
    }
    R2_lik <- 1 - exp(-2/n * (mod$logLik - ll.r))
    return(R2_lik)
}
