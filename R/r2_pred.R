#' Calculate R2_pred
#'
#' Calculate partial and total R2s for LMM, GLMM, PGLS, and PGLMM using R2_pred, an R2 based on the variance of the difference between the observed and predicted values of a fitted model.
#' 
#' @param mod A regression model with one of the following classes: 'lm', 'glm', 'lmerMod', 'glmerMod', 'phylolm', 'gls', 'pglmm', pglmm_compare', 'binaryPGLMM', or 'communityPGLMM'.
#' @param mod.r A reduced model; if not provided, the total R2 will be given by setting 'mod.r' to the model corresponding to 'mod' with the intercept as the only predictor.
#' @param phy The phylogeny for phylogenetic models (as a 'phylo' object), which must be specified for models of class `phylolm`.
#' @param gaussian.pred For models of classes `pglmm` and `pglmm_compare` when family = gaussian, which type of prediction to calculate.
#' @export
#'
#' @details  R2_pred works with classes 'lm', 'glm', 'lmerMod', 'glmerMod', 'phylolm', 'phyloglm', 'gls', 'pglmm', 'pglmm_compare', binaryPGLMM', and 'communityPGLMM' (family = gaussian and binomial).
#' 
#' \strong{LMM (lmerMod), GLMM (glmerMod), PGLMM (pglmm, pglmm_compare, binaryPGLMM and communityPGLMM):}
#' 
#' \deqn{partial R2 = 1 - var(y - y.fitted.f)/var(y - y.fitted.r)}
#' 
#' where y are the observed data, and y.fitted.f and y.fitted.r are the fitted (predicted) values from the full and reduced models. For GLMMs and PGLMMs, the values of y.fitted are in the space of the raw data (as opposed to the 'Normal' or 'latent' space). When the reduced model 'mod.r' is not specified, the total R2 is computing using the reduced model with only the intercept.
#' 
#' For pglmm and pglmm_compare with gaussian models, the default method for computing predicted values is "nearest_node" to correspond to predicted values in lmer, although the method "tip_rm" can be specified to correspond to the analyses in Ives (2018).
#' 
#' Note that the version of \code{binaryPGLMM()} in the package ape is replaced by a version contained within {rr2} that outputs all of the required information for the calculation of R2_resid.
#' 
#' \strong{PGLS (phylolm and gls):}
#' 
#' For PGLS, the total R2_pred is computed by removing each datum one at a time, predicting its value from the fitted model, repeating this for all data points, and then calculating the variance of the difference between observed and fitted values. The predictions are calculated as
#' 
#' \deqn{res.predicted[j] = V[j, -j] solve(V[-j, -j]) res[-j]}
#' 
#' where res[-j] is a vector of residuals with datum j removed, V[-j,-j] is the phylogenetic covariance matrix with row and column j removed, and V[j, -j] is row j of covariance matrix V with element j removed. The partial R2_pred is calculated from the total R2_pred from full and reduced models as
#' 
#' \deqn{partial R2 = 1 - (1 - R2_pred.f)/(1 - R2_pred.r)}
#' 
#' Note that \code{phylolm()} can have difficulties in finding solutions when there is no phylogenetic signal; when the estimate indicates no phylogenetic signal, you should refit the model with the corresponding LM.
#' 
#' \strong{LM (lm) and GLM (glm):} 
#' 
#' For compatibility and generating reduced models, rr2 will compute R2_pred for LM and GLM that correspond to LMM/PGLS and GLMM/PGLMM.
#' 
#' @author Anthony R. Ives
#' @references Ives A.R. and Li D. 2018. rr2: An R package to calculate R2s for regression models. Journal of Open Source Software. DOI:10.21105/joss.01028
#' 
#' Ives A.R. 2018. R2s for Correlated Data: Phylogenetic Models, LMMs, and GLMMs. Systematic Biology, Volume 68, Issue 2, March 2019, Pages 234-251. DOI:10.1093/sysbio/syy060
#' @seealso MuMIn, lme4, ape, phylolm, pez
#' 
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
#' R2_pred(z.f, z.x)
#' R2_pred(z.f, z.v)
#' R2_pred(z.f)
#' 
#' # GLMM with one fixed and one random effect ----
#' z.f <- glmer(y.glmm ~ x1 + (1 | u1), data = d, family = 'binomial')
#' z.x <- glmer(y.glmm ~ 1 + (1 | u1), data = d, family = 'binomial')
#' z.v <- glm(y.glmm ~ x1, data = d, family = 'binomial')
#' 
#' R2_pred(z.f, z.x)
#' R2_pred(z.f, z.v)
#' R2_pred(z.f)
#' 
#' # PGLS with a single fixed effect ----
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
#' e <- signal ^ 0.5 * rTraitCont(phy, model = 'BM', sigma = 1) +
#'   (1 - signal) ^ 0.5 * rnorm(n = n)
#' d$x <- x[match(names(e), names(x))]
#' d$y <- b1 * x + e
#' rownames(d) <- phy$tip.label
#' d$sp <- phy$tip.label
#' 
#' z.x <- gls(y ~ 1, data = d, correlation = corPagel(1, phy, form = ~sp), method = "ML")
#' z.f <- gls(y ~ x, data = d, correlation = corPagel(1, phy, form = ~sp), method = "ML")
#' z.v <- lm(y ~ x, data = d)
#' 
#' R2_pred(z.f, z.x)
#' R2_pred(z.f, z.v)
#' R2_pred(z.f)
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
#' 
#' z.f <- pglmm_compare(y ~ x, data = d, family = "binomial", phy = phy)
#' z.x <- pglmm_compare(y ~ 1, data = d, family = "binomial", phy = phy)
#' z.v <- glm(y ~ x, data = d, family = "binomial")
#' 
#' R2_pred(z.f, z.x)
#' R2_pred(z.f, z.v)
#' R2_pred(z.f)
#' 
#' #' #################
#' # A community example of pglmm {phyr} contrasting R2_pred when bayes = TRUE and bayes = F
#' 
#' library(mvtnorm)
#' nspp <- 6
#' nsite <- 4
#' 
#' # Simulate a phylogeny that has a lot of phylogenetic signal (power = 1.3)
#' phy <- compute.brlen(rtree(n = nspp), method = "Grafen", power = 1.3)
#' 
#' # Simulate species means
#' sd.sp <- 1
#' mean.sp <- rTraitCont(phy, model = "BM", sigma=sd.sp^2)
#' 
#' # Replicate values of mean.sp over sites
#' Y.sp <- rep(mean.sp, times=nsite)
#' 
#' # Simulate site means
#' sd.site <- 1
#' mean.site <- rnorm(nsite, sd=sd.site)
#' 
#' # Replicate values of mean.site over sp
#' Y.site <- rep(mean.site, each=nspp)
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
#' Y.e <- rnorm(nspp*nsite, sd = sd.e)
#' 
#' # Construct the dataset
#' d <- data.frame(sp = rep(phy$tip.label, times = nsite), site = rep(1:nsite, each = nspp))
#' 
#' # Simulate abundance data
#' d$Y <- Y.sp + Y.site + Y.attract + Y.e
#' 
#' # Full and reduced models
#' z.f <- pglmm(Y ~ 1 + (1|sp__) + (1|site) + (1|sp__@site),
#'         data = d, cov_ranef = list(sp = phy), REML = FALSE)
#' z.nested <- pglmm(Y ~ 1 + (1|sp__) + (1|site),
#'         data = d, cov_ranef = list(sp = phy), REML = FALSE)
#' z.sp <- pglmm(Y ~ 1 + (1|sp) + (1|site),
#'         data = d, cov_ranef = list(sp = phy), REML = FALSE)
#' 
#' R2_pred(z.f, z.nested)
#' R2_pred(z.nested, z.sp)
#' R2_pred(z.f)
#' 
#' # vector - matrix
#' # These are generally larger when gaussian.pred = "nearest_node"
#' R2_pred(z.f, z.nested, gaussian.pred = "nearest_node")
#' R2_pred(z.nested, z.sp, gaussian.pred = "nearest_node")
#' R2_pred(z.f, gaussian.pred = "nearest_node")
#' 
#' # # When bayes = TRUE, gaussian.pred does not work.
#' # # Commented out because INLA is not on CRAN
#' # z.f.bayes <- pglmm(Y ~ 1 + (1|sp__) + (1|site) + (1|sp__@site),
#' #                data = d, cov_ranef = list(sp = phy), bayes = TRUE)
#' # z.nested.bayes <- pglmm(Y ~ 1 + (1|sp__) + (1|site),
#' #                data = d, cov_ranef = list(sp = phy), bayes = TRUE)
#' # z.sp.bayes <- pglmm(Y ~ 1 + (1|sp) + (1|site),
#' #                data = d, cov_ranef = list(sp = phy), bayes = TRUE)
#' # 
#' # R2_pred(z.f.bayes, z.nested.bayes)
#' # R2_pred(z.nested.bayes, z.sp.bayes)
#' # R2_pred(z.f.bayes)

R2_pred <- function(mod = NULL, mod.r = NULL, gaussian.pred = "tip_rm", phy = NULL) {
    if (class(mod)[1] == "merModLmerTest") 
        class(mod) <- "lmerMod"
    
    if (!is.element(class(mod)[1], c("lm", "glm", "lmerMod", "glmerMod", "phylolm", 
                                     "gls", "pglmm", "pglmm_compare", "binaryPGLMM",
                                     "communityPGLMM"))) {
        stop("mod must be class one of classes lm, glm, lmerMod, glmerMod, phylolm (but not phyloglm), gls, pglmm, pglmm_compare, binaryPGLMM, communityPGLMM.")
    }
    
    if (class(mod)[1] == "lm") {
        if (!is.object(mod.r)) {
            y <- model.frame(mod)[, 1]
            mod.r <- lm(y ~ 1)
        }
        if (!is.element(class(mod.r)[1], c("lm"))) {
            stop("mod.r must be class lm.")
        }
        return(R2_pred.lm(mod, mod.r))
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
        return(R2_pred.glm(mod, mod.r))
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
        return(R2_pred.lmerMod(mod, mod.r))
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
        return(R2_pred.glmerMod(mod, mod.r)[1])
    }
    
    if (class(mod)[1] == "phylolm") {
        if (!is.object(mod.r)) {
            y <- mod$y
            mod.r <- lm(y ~ 1)
        }
        if (!is.element(class(mod.r)[1], c("phylolm", "lm"))) {
            stop("mod.r must be class phylolm or lm.")
        }
        return(R2_pred.phylolm(mod, mod.r, phy))
    }
    
    if (class(mod)[1] == "gls") {
      if (!is.object(mod.r)) {
        y <- as.numeric(fitted(mod)+resid(mod))
        mod.r <- lm(y ~ 1)
      }
      if (!is.element(class(mod.r)[1], c("gls", "lm"))) {
        stop("mod.r must be class gls or lm.")
      }
      return(R2_pred.gls(mod, mod.r))
    }

    if (any(class(mod) %in% c("pglmm", "communityPGLMM", "pglmm_compare"))) {
      if (!is.object(mod.r)) {
        y <- mod$Y
        if (mod$family == "gaussian") {
          mod.r <- lm(y ~ 1)
        }else{
          if(any(y < 0 | y > 1)){
            mod.r <- glm(cbind(y, mod$size-y) ~ 1, family = "binomial")
          }else{
            mod.r <- glm(y ~ 1, family = "binomial")
          }
        }
      }
      
      if (!any(class(mod.r) %in% c("pglmm", "communityPGLMM", "pglmm_compare", "lm", "glm"))) {
        stop("mod.r must be of class pglmm, communityPGLMM, pglmm_compare, lm, or glm.")
      }
      if (mod$family == "gaussian") {
        if(gaussian.pred == "nearest_node") warning("Predictions are made with gaussian.pred = nearest_node")
        if(mod$bayes == TRUE) warning("Predictions are made with gaussian.pred = nearest_node for models with bayes = TRUE")
        return(R2_pred.pglmm.gaussian(mod, mod.r, gaussian.pred = gaussian.pred))
      }
        
      if (mod$family != "gaussian") 
        return(R2_pred.pglmm.glm(mod, mod.r))
    }

    if (class(mod)[1] == "binaryPGLMM") {
      if (!is.object(mod.r)) {
        y <- mod$y
        if(any(y < 0 | y > 1)){
          mod.r <- glm(cbind(y, mod$size-y) ~ 1, family = "binomial")
        }else{
          mod.r <- glm(y ~ 1, family = "binomial")
        }
      }
      if (!is.element(class(mod.r)[1], c("binaryPGLMM", "glm"))) {
        stop("mod.r must be class binaryPGLMM or glm.")
      }
      return(R2_pred.binaryPGLMM(mod, mod.r)[1])
    }
    
    # if (class(mod)[1] == "communityPGLMM") {
    #     if (!is.object(mod.r)) {
    #         y <- mod$Y
    #         if (mod$family == "gaussian") {
    #             mod.r <- lm(y ~ 1)
    #         }
    #         if (mod$family == "binomial") {
    #             mod.r <- glm(y ~ 1, family = "binomial")
    #         }
    #     }
    #     
    #     if (any(class(mod) %in% c("communityPGLMM", "lm", "glm"))) {
    #         stop("mod.r must be of class communityPGLMM, lm, or glm (binomial model).")
    #     }
    #     
    #     if (mod$family == "gaussian") 
    #         return(R2_pred.communityPGLMM.gaussian(mod, mod.r))
    #     if (mod$family == "binomial") 
    #         return(R2_pred.communityPGLMM.binomial(mod, mod.r))
    # }
}

R2_pred.lm <- function(mod = NA, mod.r = NA) {
    y <- model.frame(mod)[, 1]
    SSE.pred <- var(y - stats::fitted(mod))
    SSE.pred.r <- var(y - stats::fitted(mod.r))
    R2_pred <- 1 - SSE.pred/SSE.pred.r
    return(R2_pred)
}

R2_pred.glm <- function(mod = NA, mod.r = NA) {
  y <- mod$y
  SSE.pred <- var(y - stats::fitted(mod))
  SSE.pred.r <- var(y - stats::fitted(mod.r))
  R2_pred <- 1 - SSE.pred/SSE.pred.r
  return(R2_pred)
}

R2_pred.lmerMod <- function(mod = NA, mod.r = NA) {
    y <- model.frame(mod)[, 1]
    SSE.pred <- var(y - fitted(mod))
    SSE.pred.r <- var(y - stats::fitted(mod.r))
    R2_pred <- 1 - SSE.pred/SSE.pred.r
    return(R2_pred)
}

R2_pred.glmerMod <- function(mod = NA, mod.r = NA) {
  y <- model.frame(mod)[, 1]
  if(is.matrix(y)) y <- y[,1]/rowSums(y)
  SSE.pred <- var(y - fitted(mod))
  SSE.pred.r <- var(y - stats::fitted(mod.r))
  R2_pred <- 1 - SSE.pred/SSE.pred.r
  return(R2_pred)
}

R2_pred.phylolm <- function(mod = NULL, mod.r = NULL, phy = NULL) {
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
  if (inherits(mod.r, "phylolm")) {
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
  
  if (inherits(mod.r, "lm")) {
    Yhat.r <- stats::fitted(mod.r)
  }
  
  SSE.pred.r <- var(y - Yhat.r)
  
  return(1 - SSE.pred/SSE.pred.r)
}

R2_pred.gls <- function(mod = NULL, mod.r = NULL) {
  y <- as.numeric(fitted(mod.r)+resid(mod.r))
  n <- mod$dims$N
  
  V <- nlme::corMatrix(mod$modelStruct$corStruct)
  if(length(V)>1) {
    V <- Matrix::bdiag(V)
    V <- as.matrix(V)
  }
  if(!is.null(attr(mod$modelStruct$varStruct, 'weights'))){
    Vdiag <- 1/attr(mod$modelStruct$varStruct, 'weights')
    V <- diag(Vdiag) %*% V %*% diag(Vdiag)
  }
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
  if (inherits(mod.r, "gls")) {

    V.r <- nlme::corMatrix(mod.r$modelStruct$corStruct)
    if(length(V.r)>1) {
      V.r <- Matrix::bdiag(V.r)
      V.r <- as.matrix(V.r)
    }
    if(!is.null(attr(mod.r$modelStruct$varStruct, 'weights'))){
      Vdiag.r <- 1/attr(mod.r$modelStruct$varStruct, 'weights')
      V.r <- diag(Vdiag.r) %*% V.r %*% diag(Vdiag.r)
    }
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
  
  if (inherits(mod.r, "lm")) {
    Yhat.r <- stats::fitted(mod.r)
  }
  
  SSE.pred.r <- var(y - Yhat.r)
  
  return(1 - SSE.pred/SSE.pred.r)
}

R2_pred.binaryPGLMM <- function(mod = NULL, mod.r = NULL) {
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

# pglmm.gaussian.predict <- function(mod) {
#     y <- matrix(mod$Y, ncol=1)
#     X <- mod$X
#     n <- dim(X)[1]
#     fit <- X %*% mod$B
#     R <- y - fit
#     V <- solve(mod$iV)
#     v <- V
#     for (i in 1:n) {
#         v[i, i] <- max(V[i, -i])
#     }
#     Rhat <- v %*% mod$iV %*% R
#     Yhat <- as.numeric(fit + Rhat)
#     return(Yhat)
# }

R2_pred.pglmm.glm <- function(mod = NULL, mod.r = NULL) {
  # full model
  SSE.pred <- var(mod$Y - mod$mu)

  # reduced model
  if (any(is.element(class(mod.r), c("pglmm", "pglmm_compare","communityPGLMM")))){
    SSE.pred.r <- var(mod.r$Y - mod.r$mu)
  }
  if (inherits(mod.r, "glm")) {
    Yhat.r <- stats::fitted(mod.r)
    SSE.pred.r <- var(mod$Y - Yhat.r)
  }
  
  return(as.numeric(1 - SSE.pred/SSE.pred.r))
}

R2_pred.pglmm.gaussian <- function(mod = NULL, mod.r = NULL, gaussian.pred = gaussian.pred) {

    # full model
    Yhat <- pglmm_predicted_values(mod, gaussian.pred = gaussian.pred)
    SSE.pred <- as.numeric(var(mod$Y - Yhat))

    # reduced model
    if (any(is.element(class(mod.r), c("pglmm", "pglmm_compare","communityPGLMM")))){
        Yhat.r <- pglmm_predicted_values(mod.r, gaussian.pred = gaussian.pred)
        SSE.pred.r <- as.numeric(var(mod.r$Y - Yhat.r))
    }
    
    if (inherits(mod.r, "lm")) {
        y.r <- model.frame(mod.r)[, 1]
        Yhat.r <- stats::fitted(mod.r)
        SSE.pred.r <- var(y.r - Yhat.r)
    }
    
    return(as.numeric(1 - SSE.pred/SSE.pred.r))
}

# R2_pred.communityPGLMM.binomial <- function(mod = NULL, mod.r = NULL) {
#     y <- mod$Y
#     Yhat <- mod$mu
#     SSE.pred <- var(y - Yhat)
#     
#     if (class(mod.r)[1] == "communityPGLMM") {
#         y.r <- mod.r$Y
#         Yhat.r <- mod.r$mu
#         SSE.pred.r <- var(y.r - Yhat.r)
#     }
#     
#     if (class(mod.r)[1] == "glm") {
#         y.r <- mod.r$y
#         Yhat.r <- mod.r$fitted.values
#         SSE.pred.r <- var(y.r - Yhat.r)
#     }
#     
#     return(1 - SSE.pred/SSE.pred.r)
# }
# 
# # these two versions for communitypglmm and lmer work like R2_pred for phylolm
# # objects, in which the points are predicted for tips after removing the tip
# # values.
# pglmm_predict_alt <- function(mod) {
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
# R2_pred.communityPGLMM.gaussian.alt <- function(mod = NULL, mod.r = NULL) {
#     Yhat <- pglmm_predict_alt(mod)
#     SSE.pred <- var(mod$Y - Yhat)
#     
#     # reduced model
#     if (class(mod.r) == "communityPGLMM") {
#         Yhat.r <- pglmm_predict_alt(mod.r)
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
# R2_pred.lmer.gaussian.alt <- function(mod = NULL, mod.r = NULL) {
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
