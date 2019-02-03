#' Calculate R2.pred
#'
#' Calculate partial and total R2s for LMM, GLMM, PGLS, and PGLMM using R2.pred, an R2 based on the variance of the difference between the observed and predicted values of a fitted model.
#' 
#' @param mod A regression model with one of the following classes: 'lm', 'glm', 'lmerMod', 'glmerMod', 'phylolm', 'gls', 'binaryPGLMM', or 'communityPGLMM'.
#' @param mod.r A reduced model; if not provided, the total R2 will be given by setting 'mod.r' to the model corresponding to 'mod' with the intercept as the only predictor.
#' @param phy The phylogeny for phylogenetic models (as a 'phylo' object), which must be specified for models of class `phylolm`.
#' @return R2.pred value.
#' @export
#'
#' @details  R2.pred works with classes 'lm', 'glm', 'lmerMod', 'glmerMod', 'phylolm', 'phyloglm', 'gls', binaryPGLMM', and 'communityPGLMM' (family = gaussian and binomial).
#' 
#' \strong{LMM (lmerMod), GLMM (glmerMod), PGLMM (binaryPGLMM and communityPGLMM):}
#' 
#' \deqn{partial R2 = 1 - var(y - y.fitted.f)/var(y - y.fitted.r)}
#' 
#' where y are the observed data, and y.fitted.f and y.fitted.r are the fitted (predicted) values from the full and reduced models. For GLMMs and PGLMMs, the values of y.fitted are in the space of the raw data (as opposed to the 'Normal' or 'latent' space). When the reduced model 'mod.r' is not specified, the total R2 is computing using the reduced model with only the intercept.
#' 
#' Note that the version of \code{binaryPGLMM()} in the package ape is replaced by a version contained within {rr2} that outputs all of the required information for the calculation of R2.resid.
#' 
#' \strong{PGLS (phylolm and gls):}
#' 
#' For PGLS, the total R2.pred is computed by removing each datum one at a time, predicting its value from the fitted model, repeating this for all data points, and then calculating the variance of the difference between observed and fitted values. The predictions are calculated as
#' 
#' \deqn{res.predicted[j] = V[j, -j] solve(V[-j, -j]) res[-j]}
#' 
#' where res[-j] is a vector of residuals with datum j removed, V[-j,-j] is the phylogenetic covariance matrix with row and column j removed, and V[j, -j] is row j of covariance matrix V with element j removed. The partial R2.pred is calculated from the total R2.pred from full and reduced models as
#' 
#' \deqn{partial R2 = 1 - (1 - R2.pred.f)/(1 - R2.pred.r)}
#' 
#' Note that \code{phylolm()} can have difficulties in finding solutions when there is no phylogenetic
#' signal; when the estimate indicates no phylogenetic signal, you should refit the model with the corresponding LM.
#' 
#' \strong{LM (lm) and GLM (glm):} 
#' 
#' For compatibility and generating reduced models, rr2 will compute R2.pred for LM and GLM that correspond to LMM/PGLS and GLMM/PGLMM.
#' 
#' @author Anthony R. Ives
#' @references Ives A.R. and Li D. 2018. rr2: An R package to calculate R2s for regression models. Journal of Open Source Software. DOI:10.21105/joss.01028
#' 
#' Ives A.R. 2018. R2s for Correlated Data: Phylogenetic Models, LMMs, and GLMMs. Systematic Biology. DOI:10.1093/sysbio/syy060
#' @seealso MuMIn, lme4, ape, phylolm, pez
#' 
#' @examples library(ape)
#' library(phylolm)
#' library(lme4)
#' library(nlme)
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
#' R2.pred(z.f, z.x)
#' R2.pred(z.f, z.v)
#' R2.pred(z.f)
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
#' R2.pred(z.f, z.x)
#' R2.pred(z.f, z.v)
#' R2.pred(z.f)
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
#' z.x <- gls(y ~ 1, data = d, correlation = corPagel(1, phy), method = "ML")
#' z.f <- gls(y ~ x, data = d, correlation = corPagel(1, phy), method = "ML")
#' z.v <- lm(y ~ x, data = d)
#' 
#' R2.pred(z.f, z.x)
#' R2.pred(z.f, z.v)
#' R2.pred(z.f)
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
#' # Use the function binaryPGLMM() from the rr2 package rather than ape.
#' z.f <- rr2::binaryPGLMM(y ~ x, data = d, phy = phy)
#' z.x <- rr2::binaryPGLMM(y ~ 1, data = d, phy = phy)
#' z.v <- glm(y ~ x, data = d, family = 'binomial')
#' 
#' R2.pred(z.f, z.x)
#' R2.pred(z.f, z.v)
#' R2.pred(z.f)
#' 
R2.pred <- function(mod = NULL, mod.r = NULL, phy = NULL) {
    if (class(mod)[1] == "merModLmerTest") 
        class(mod) <- "lmerMod"
    
    if (!is.element(class(mod)[1], c("lm", "glm", "lmerMod", "glmerMod", "phylolm", "gls", "binaryPGLMM", "communityPGLMM"))) {
        stop("mod must be class one of classes lm, glm, lmerMod, glmerMod, phylolm (but not phyloglm), gls, binaryPGLMM, communityPGLMM.")
    }
    
    if (class(mod)[1] == "lm") {
        if (!is.object(mod.r)) {
            y <- model.frame(mod)[, 1]
            mod.r <- lm(y ~ 1)
        }
        if (!is.element(class(mod.r)[1], c("lm"))) {
            stop("mod.r must be class lm.")
        }
        return(R2.pred.lm(mod, mod.r))
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
        return(R2.pred.glm(mod, mod.r))
    }
    
    if (class(mod)[1] == "lmerMod") {
        if (!is.object(mod.r)) {
            y <- model.frame(mod)[, 1]
            mod.r <- lm(y ~ 1)
        }
        if (class(mod.r)[1] == "merModLmerTest") 
            class(mod.r) <- "lmerMod"
        if (!is.element(class(mod)[1], c("lmerMod", "lm"))) {
            stop("mod.r must be class lmerMod or lm.")
        }
        return(R2.pred.lmerMod(mod, mod.r))
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
        return(R2.pred.glmerMod(mod, mod.r)[1])
    }
    
    if (class(mod)[1] == "phylolm") {
        if (!is.object(phy)) {
            stop("For phylolm you must provide the phylo object.")
        }
        if (!is.object(mod.r)) {
            y <- mod$y
            mod.r <- lm(y ~ 1)
        }
        if (!is.element(class(mod.r)[1], c("phylolm", "lm"))) {
            stop("mod.r must be class phylolm or lm.")
        }
        return(R2.pred.phylolm(mod, mod.r, phy))
    }
    
    if (class(mod)[1] == "gls") {
      if (!is.object(mod.r)) {
        y <- as.numeric(fitted(mod)+resid(mod))
        mod.r <- lm(y ~ 1)
      }
      if (!is.element(class(mod.r)[1], c("gls", "lm"))) {
        stop("mod.r must be class gls or lm.")
      }
      return(R2.pred.gls(mod, mod.r))
    }
    
    if (class(mod)[1] == "binaryPGLMM") {
        if (!is.object(mod.r)) {
            y <- mod$y
            mod.r <- glm(y ~ 1, family = "binomial")
        }
        if (!is.element(class(mod.r)[1], c("binaryPGLMM", "glm"))) {
            stop("mod.r must be class binaryPGLMM or glm.")
        }
        return(R2.pred.binaryPGLMM(mod, mod.r)[1])
    }
    
    if (class(mod)[1] == "communityPGLMM") {
        if (!is.object(mod.r)) {
            y <- mod$Y
            if (mod$family == "gaussian") {
                mod.r <- lm(y ~ 1)
            }
            if (mod$family == "binomial") {
                mod.r <- glm(y ~ 1, family = "binomial")
            }
        }
        
        if (!is.element(class(mod.r)[1], c("communityPGLMM", "lm", "glm"))) {
            stop("mod.r must be of class communityPGLMM, lm, or glm (binomial model).")
        }
        
        if (mod$family == "gaussian") 
            return(R2.pred.communityPGLMM.gaussian(mod, mod.r))
        if (mod$family == "binomial") 
            return(R2.pred.communityPGLMM.binomial(mod, mod.r))
    }
}

R2.pred.lm <- function(mod = NA, mod.r = NA) {
    y <- model.frame(mod)[, 1]
    SSE.pred <- var(y - stats::fitted(mod))
    SSE.pred.r <- var(y - stats::fitted(mod.r))
    R2.pred <- 1 - SSE.pred/SSE.pred.r
    return(R2.pred)
}

R2.pred.glm <- function(mod = NA, mod.r = NA) {
  y <- model.frame(mod)[, 1]
  if(is.matrix(y)) y <- y[,1]/rowSums(y)
  SSE.pred <- var(y - stats::fitted(mod))
  SSE.pred.r <- var(y - stats::fitted(mod.r))
  R2.pred <- 1 - SSE.pred/SSE.pred.r
  return(R2.pred)
}

R2.pred.lmerMod <- function(mod = NA, mod.r = NA) {
    y <- model.frame(mod)[, 1]
    SSE.pred <- var(y - fitted(mod))
    SSE.pred.r <- var(y - stats::fitted(mod.r))
    R2.pred <- 1 - SSE.pred/SSE.pred.r
    return(R2.pred)
}

R2.pred.glmerMod <- function(mod = NA, mod.r = NA) {
  y <- model.frame(mod)[, 1]
  if(is.matrix(y)) y <- y[,1]/rowSums(y)
  SSE.pred <- var(y - fitted(mod))
  SSE.pred.r <- var(y - stats::fitted(mod.r))
  R2.pred <- 1 - SSE.pred/SSE.pred.r
  return(R2.pred)
}

R2.pred.phylolm <- function(mod = NULL, mod.r = NULL, phy = NULL) {
  y <- mod$y
  X <- mod$X
  n <- dim(X)[1]
  
  if (!mod$model %in% c("lambda", "OUrandomRoot", "OUfixedRoot", "BM", "kappa", 
                        "delta", "EB", "trend")) {
    stop("Evolution model not supported yet.")
  }
  
  phy.f <- transf_phy(mod, phy)  # function in the utils.R
  
  V <- ape::vcv(phy.f)
  R <- y - stats::fitted(mod)
  
  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    VV <- V[-j, -j]
    iVV <- solve(VV)
    v <- V[j, -j]
    Rhat[j] <- v %*% iVV %*% r
  }
  
  Yhat <- as.numeric(stats::fitted(mod) + Rhat)
  SSE.pred <- var(y - Yhat)
  
  # reduced model
  if (class(mod.r) == "phylolm") {
    if (!mod.r$model %in% c("lambda", "OUrandomRoot", "OUfixedRoot", "BM", "kappa", 
                            "delta", "EB", "trend")) {
      stop("Evolution model not supported yet.")
    }
    
    phy.r <- transf_phy(mod.r, phy)
    
    V.r <- ape::vcv(phy.r)
    R.r <- y - stats::fitted(mod.r)
    Rhat.r <- matrix(0, nrow = n, ncol = 1)
    for (j in 1:n) {
      r.r <- R.r[-j]
      VV.r <- V.r[-j, -j]
      iVV.r <- solve(VV.r)
      v.r <- V.r[j, -j]
      Rhat.r[j] <- v.r %*% iVV.r %*% r.r
    }
    Yhat.r <- as.numeric(stats::fitted(mod.r) + Rhat.r)
  }
  
  if (class(mod.r) == "lm") {
    Yhat.r <- stats::fitted(mod.r)
  }
  
  SSE.pred.r <- var(y - Yhat.r)
  
  return(1 - SSE.pred/SSE.pred.r)
}

R2.pred.gls <- function(mod = NULL, mod.r = NULL) {
  y <- as.numeric(fitted(mod.r)+resid(mod.r))
  n <- mod$dims$N
  
  V <- nlme::corMatrix(mod$modelStruct$corStruct)
  R <- y - stats::fitted(mod)
  
  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    VV <- V[-j, -j]
    iVV <- solve(VV)
    v <- V[j, -j]
    Rhat[j] <- v %*% iVV %*% r
  }
  
  Yhat <- as.numeric(fitted(mod) + Rhat)
  SSE.pred <- var(y - Yhat)
  
  # reduced model
  if (class(mod.r) == "gls") {

    V.r <- nlme::corMatrix(mod.r$modelStruct$corStruct)
    R.r <- y - fitted(mod.r)
    Rhat.r <- matrix(0, nrow = n, ncol = 1)
    for (j in 1:n) {
      r.r <- R.r[-j]
      VV.r <- V.r[-j, -j]
      iVV.r <- solve(VV.r)
      v.r <- V.r[j, -j]
      Rhat.r[j] <- v.r %*% iVV.r %*% r.r
    }
    Yhat.r <- as.numeric(fitted(mod.r) + Rhat.r)
  }
  
  if (class(mod.r) == "lm") {
    Yhat.r <- stats::fitted(mod.r)
  }
  
  SSE.pred.r <- var(y - Yhat.r)
  
  return(1 - SSE.pred/SSE.pred.r)
}

R2.pred.binaryPGLMM <- function(mod = NULL, mod.r = NULL) {
    Yhat <- mod$mu
    y <- mod$y
    SSE.pred <- var(y - Yhat)
    
    if (class(mod.r)[1] == "binaryPGLMM") {
        Yhat.r <- mod.r$mu
        SSE.pred.r <- var(y - Yhat.r)
    }
    
    if (class(mod.r)[1] == "glm") {
        Yhat.r <- mod.r$fitted.values
        SSE.pred.r <- var(y - Yhat.r)
    }
    
    return(1 - SSE.pred/SSE.pred.r)
}

pglmm.predict <- function(mod) {
    y <- mod$Y
    X <- mod$X
    n <- dim(X)[1]
    fit <- X %*% mod$B
    R <- y - fit
    V <- solve(mod$iV)
    v <- y
    for (i in 1:n) {
        v[i, i] <- max(V[i, -i])
    }
    Rhat <- v %*% mod$iV %*% R
    Yhat <- as.numeric(fit + Rhat)
    return(Yhat)
}

R2.pred.communityPGLMM.gaussian <- function(mod = NULL, mod.r = NULL) {
    Yhat <- pglmm.predict(mod)
    SSE.pred <- var(mod$Y - Yhat)
    
    # reduced model
    if (class(mod.r) == "communityPGLMM") {
        Yhat.r <- pglmm.predict(mod.r)
        SSE.pred.r <- var(mod.r$Y - Yhat.r)
    }
    
    if (class(mod.r) == "lm") {
        y.r <- model.frame(mod.r)[, 1]
        Yhat.r <- stats::fitted(mod.r)
        SSE.pred.r <- var(y.r - Yhat.r)
    }
    
    return(1 - SSE.pred/SSE.pred.r)
}

R2.pred.communityPGLMM.binomial <- function(mod = NULL, mod.r = NULL) {
    y <- mod$Y
    Yhat <- mod$mu
    SSE.pred <- var(y - Yhat)
    
    if (class(mod.r)[1] == "communityPGLMM") {
        y.r <- mod.r$Y
        Yhat.r <- mod.r$mu
        SSE.pred.r <- var(y.r - Yhat.r)
    }
    
    if (class(mod.r)[1] == "glm") {
        y.r <- mod.r$y
        Yhat.r <- mod.r$fitted.values
        SSE.pred.r <- var(y.r - Yhat.r)
    }
    
    return(1 - SSE.pred/SSE.pred.r)
}
# 
# # these two versions for communitypglmm and lmer work like r2.pred for phylolm
# # objects, in which the points are predicted for tips after removing the tip
# # values.
# pglmm.predict.alt <- function(mod) {
#     Y <- mod$Y
#     X <- mod$X
#     n <- dim(X)[1]
#     fit <- X %*% mod$B
#     R <- Y - fit
#     V <- solve(mod$iV)
#     Rhat <- matrix(0, nrow = n, ncol = 1)
#     for (j in 1:n) {
#         r <- R[-j]
#         VV <- V[-j, -j]
#         iVV <- solve(VV)
#         v <- V[j, -j]
#         Rhat[j] <- v %*% iVV %*% r
#     }
#     Yhat <- as.numeric(fit + Rhat)
#     return(Yhat)
# }
# 
# R2.pred.communityPGLMM.gaussian.alt <- function(mod = NULL, mod.r = NULL) {
#     Yhat <- pglmm.predict.alt(mod)
#     SSE.pred <- var(mod$Y - Yhat)
#     
#     # reduced model
#     if (class(mod.r) == "communityPGLMM") {
#         Yhat.r <- pglmm.predict.alt(mod.r)
#         SSE.pred.r <- var(mod.r$Y - Yhat.r)
#     }
#     
#     if (class(mod.r) == "lm") {
#         Y.r <- model.frame(mod.r)[, 1]
#         Yhat.r <- stats::fitted(mod.r)
#         SSE.pred.r <- var(Y.r - Yhat.r)
#     }
#     
#     return(1 - SSE.pred/SSE.pred.r)
# }
# 
# mer.predict.alt <- function(mod) {
#     Y <- model.frame(mod)[, 1]
#     X <- model.matrix(mod)
#     n <- dim(X)[1]
#     fit <- X %*% lme4::fixef(mod)
#     R <- Y - fit
#     s2 <- lme4::getME(mod, "theta")^2
#     n.re <- length(s2)
#     Zt <- lme4::getME(mod, "Ztlist")
#     V <- diag(n)
#     for (i in 1:n.re) V <- V + s2[i] * t(Zt[[i]]) %*% Zt[[i]]
#     V <- sigma(mod)^2 * V
#     V <- V/det(V)^(1/n)
#     
#     Rhat <- matrix(0, nrow = n, ncol = 1)
#     for (j in 1:n) {
#         r <- R[-j]
#         VV <- V[-j, -j]
#         iVV <- solve(VV)
#         v <- V[j, -j]
#         Rhat[j] <- v %*% iVV %*% r
#     }
#     
#     Yhat <- as.numeric(fit + Rhat)
# }
# 
# R2.pred.lmer.gaussian.alt <- function(mod = NULL, mod.r = NULL) {
#     Y <- model.frame(mod)[, 1]
#     Yhat <- mer.predict.alt(mod)
#     SSE.pred <- var(Y - Yhat)
#     
#     # reduced model
#     if (class(mod.r) == "lmerMod") {
#         Yhat.r <- mer.predict.alt(mod.r)
#     }
#     
#     if (class(mod.r) == "lm") {
#         Yhat.r <- stats::fitted(mod.r)
#     }
#     
#     SSE.pred.r <- var(Y - Yhat.r)
#     return(1 - SSE.pred/SSE.pred.r)
# }
