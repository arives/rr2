#' Calculate R2.pred
#'
#' Calculate R2.pred for LMM, GLMM, PGLM, and PGLMMs.
#' @param mod a regression model with the following class: 'lmerMod', 'glmerMod', 'phylolm', and 'binaryPGLMM'
#' @param mod.r reduced model, if not provided, will use corresponding models with intercept as the only predictor
#' @param phy a phylogeny with tip labels and branch length
#' @return R2.pred
#' @export
#'
R2.pred <- function(mod = NULL, mod.r = NULL, phy = NULL) {
  
  if (!is.element(class(mod)[1], c("lm", "glm", "lmerMod", "glmerMod", "phylolm", "binaryPGLMM", "communityPGLMM"))) {
    stop("mod must be class one of classes lm, glm, lmerMod, glmerMod, phylolm (but not phyloglm), binaryPGLMM, communityPGLMM.")
  }
  
  if (class(mod)[1] == "lm") {
    if (!is.object(mod.r)) {
      Y <- model.frame(mod)[, 1]
      mod.r <- lm(Y ~ 1)
    }
    if (!is.element(class(mod.r)[1], c("lm"))) {
      stop("mod.r must be class lm.")
    }
    return(R2.pred.lm(mod, mod.r))
  }
  
  if (class(mod)[1] == "glm") {
    if (!is.object(mod.r)) {
      Y <- model.frame(mod)[, 1]
      mod.r <- glm(Y ~ 1, family = family(mod)[[1]])
    }
    if (!is.element(class(mod.r)[1], c("glm"))) {
      stop("mod.r must be class glm.")
    }
    return(R2.pred.glm(mod, mod.r))
  }
  
  if (class(mod)[1] == "lmerMod") {
    if (!is.object(mod.r)) {
      Y <- model.frame(mod)[, 1]
      mod.r <- lm(Y ~ 1)
    }
    if (!is.element(class(mod)[1], c("lmerMod", "lm"))) {
      stop("mod.r must be class lmerMod or lm.")
    }
    return(R2.pred.lmerMod(mod, mod.r))
  }
  
  if (class(mod)[1] == "glmerMod") {
    if (family(mod)[[1]] != "binomial") {
      stop("Sorry, but only binomial (binary) glmerMod models are allowed.")
    }
    if (!is.object(mod.r)) {
      Y <- model.frame(mod)[, 1]
      mod.r <- glm(Y ~ 1, family = family(mod)[[1]])
    }
    if (!is.element(class(mod.r)[1], c("glmerMod", "glm"))) {
      stop("mod.r must be class glmerMod or glm.")
    }
    return(R2.pred.glmerMod(mod, mod.r)[1])
  }
  
  if (class(mod)[1] == "phylolm") {
    if (!is.object(phy)) {
      stop("For phylolm you must provide the phylo object")
    }
    if (!is.object(mod.r)) {
      Y <- mod$y
      mod.r <- lm(Y ~ 1)
    }
    if (!is.element(class(mod.r)[1], c("phylolm", "lm"))) {
      stop("mod.r must be class phylolm or lm.")
    }
    return(R2.pred.phylolm(mod, mod.r, phy))
  }
  
  if (class(mod)[1] == "binaryPGLMM") {
    if (!is.object(mod.r)) {
      Y <- mod$y
      mod.r <- glm(Y ~ 1, family = "binomial")
    }
    if (!is.element(class(mod.r)[1], c("binaryPGLMM", "glm"))) {
      stop("mod.r must be class binaryPGLMM or glm.")
    }
    return(R2.pred.binaryPGLMM(mod, mod.r)[1])
  }
  
  if (class(mod)[1] == "communityPGLMM") {
    if(!is.object(mod.r)){
      Y <- mod$Y
      if(mod$family == "gaussian"){
        mod.r <- lm(Y ~ 1)
      }
      if(mod$family == "binomial"){
        mod.r <- glm(Y ~ 1, family = "binomial")
      }
    }
    
    if (!is.element(class(mod.r)[1], c("communityPGLMM", "lm", "glm"))) {
      stop("mod.r must be of class communityPGLMM, lm, or glm (binomial model).")
    }
    
    if(mod$family == "gaussian") return(R2.pred.communityPGLMM.gaussian(mod, mod.r))
    if(mod$family == "binomial") return(R2.pred.communityPGLMM.binomial(mod, mod.r))
  }
}

R2.pred.lm <- function(mod = NA, mod.r = NA) {
  Y <- model.frame(mod)[, 1]
  SSE.pred <- var(Y - stats::fitted(mod))
  SSE.pred.r <- var(Y - stats::fitted(mod.r))
  R2.pred <- 1 - SSE.pred/SSE.pred.r
  return(R2.pred)
}

R2.pred.glm <- R2.pred.lm

R2.pred.lmerMod <- function(mod = NA, mod.r = NA) {
  Y <- model.frame(mod)[,1]
  SSE.pred <- var(Y - fitted(mod))
  SSE.pred.r <- var(Y - stats::fitted(mod.r))
  R2.pred <- 1 - SSE.pred/SSE.pred.r
  return(R2.pred)
}

R2.pred.glmerMod <- R2.pred.lmerMod

R2.pred.phylolm <- function(mod = NULL, mod.r = NULL, phy = NULL) {
  Y <- mod$y
  X <- mod$X
  n <- dim(X)[1]
  
  if (!mod$model %in% c("lambda", "OUrandomRoot", "OUfixedRoot", "BM", "kappa", "delta", "EB", "trend")) {
    stop("evolution model not supported yet")
  }
  
  phy.f <- transf_phy(mod, phy) # function in the utils.R
  
  V <- ape::vcv(phy.f)
  R <- Y - stats::fitted(mod)
  
  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    VV <- V[-j, -j]
    iVV <- solve(VV)
    v <- V[j, -j]
    Rhat[j] <- v %*% iVV %*% r
  }
  
  Yhat <- as.numeric(stats::fitted(mod) + Rhat)
  SSE.pred <- var(Y - Yhat)
  
  # reduced model
  if (class(mod.r) == "phylolm") {
    if (!mod.r$model %in% c("lambda", "OUrandomRoot", "OUfixedRoot", "BM", "kappa", "delta", "EB", "trend")) {
      stop("evolution model not supported yet")
    }
    
    phy.r <- transf_phy(mod.r, phy)
    
    V.r <- ape::vcv(phy.r)
    R.r <- Y - stats::fitted(mod.r)
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
  
  SSE.pred.r <- var(Y - Yhat.r)
  
  return(1 - SSE.pred/SSE.pred.r)
}

R2.pred.binaryPGLMM <- function(mod = NULL, mod.r = NULL) {
  Yhat <- mod$mu
  Y <- mod$y
  SSE.pred <- var(Y - Yhat)
  
  if (class(mod.r)[1] == "binaryPGLMM") {
    Yhat.r <- mod.r$mu
    SSE.pred.r <- var(Y - Yhat.r)
  }
  
  if (class(mod.r)[1] == "glm") {
    Yhat.r <- mod.r$fitted.values
    SSE.pred.r <- var(Y - Yhat.r)
  }
  
  return(1 - SSE.pred/SSE.pred.r)
}

pglmm.predict <- function(mod){
  Y <- mod$Y
  X <- mod$X
  n <- dim(X)[1]
  fit <- X %*% mod$B
  R <- Y - fit
  V <- solve(mod$iV)
  v <- V
  for(i in 1:n) {
    v[i,i] <- max(V[i, -i])
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
    Y.r <- model.frame(mod.r)[, 1]
    Yhat.r <- stats::fitted(mod.r)
    SSE.pred.r <- var(Y.r - Yhat.r)
  }
  
  return(1 - SSE.pred/SSE.pred.r)
}

R2.pred.communityPGLMM.binomial <- function(mod = NULL, mod.r = NULL) {
  Y <- mod$Y
  Yhat <- mod$mu
  SSE.pred <- var(Y - Yhat)
  
  if (class(mod.r)[1] == "communityPGLMM") {
    Y.r <- mod.r$Y
    Yhat.r <- mod.r$mu
    SSE.pred.r <- var(Y.r - Yhat.r)
  }
  
  if(class(mod.r)[1] == "glm"){
    Y.r <- mod.r$y
    Yhat.r <- mod.r$fitted.values
    SSE.pred.r <- var(Y.r - Yhat.r)
  }
  
  return(1 - SSE.pred/SSE.pred.r)
}

# these two versions for communitypglmm and lmer work like r2.pred for phylolm objects, 
# in which the points are predicted for tips after removing the tip values.
pglmm.predict.alt <- function(mod){
  Y <- mod$Y
  X <- mod$X
  n <- dim(X)[1]
  fit <- X %*% mod$B
  R <- Y - fit
  V <- solve(mod$iV)
  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    VV <- V[-j, -j]
    iVV <- solve(VV)
    v <- V[j, -j]
    Rhat[j] <- v %*% iVV %*% r
  }
  Yhat <- as.numeric(fit + Rhat)
  return(Yhat)
}

R2.pred.communityPGLMM.gaussian.alt <- function(mod = NULL, mod.r = NULL) {
  Yhat <- pglmm.predict.alt(mod)
  SSE.pred <- var(mod$Y - Yhat)
  
  # reduced model
  if (class(mod.r) == "communityPGLMM") {
    Yhat.r <- pglmm.predict.alt(mod.r)
    SSE.pred.r <- var(mod.r$Y - Yhat.r)
  }
  
  if (class(mod.r) == "lm") {
    Y.r <- model.frame(mod.r)[, 1]
    Yhat.r <- stats::fitted(mod.r)
    SSE.pred.r <- var(Y.r - Yhat.r)
  }

  return(1 - SSE.pred/SSE.pred.r)
}

mer.predict.alt <- function(mod){
  Y <- model.frame(mod)[,1]
  X <- model.matrix(mod)
  n <- dim(X)[1]
  fit <- X %*% lme4::fixef(mod)
  R <- Y - fit
  s2 <- lme4::getME(mod,"theta")^2
  n.re <- length(s2)
  Zt <- lme4::getME(mod, "Ztlist")
  V <- diag(n)
  for(i in 1:n.re) V <- V + s2[i] * t(Zt[[i]]) %*% Zt[[i]]
  V <- sigma(mod)^2 * V
  V <- V/det(V)^(1/n)
  
  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    VV <- V[-j, -j]
    iVV <- solve(VV)
    v <- V[j, -j]
    Rhat[j] <- v %*% iVV %*% r
  }
  
  Yhat <- as.numeric(fit + Rhat)
}

R2.pred.lmer.gaussian.alt <- function(mod = NULL, mod.r = NULL) {
  Y <- model.frame(mod)[, 1]
  Yhat <- mer.predict.alt(mod)
  SSE.pred <- var(Y - Yhat)
  
  # reduced model
  if (class(mod.r) == "lmerMod") {
    Yhat.r <- mer.predict.alt(mod.r)
  }
  
  if (class(mod.r) == "lm") {
    Yhat.r <- stats::fitted(mod.r)
  }
  
  SSE.pred.r <- var(Y - Yhat.r)
  return(1 - SSE.pred/SSE.pred.r)
}
