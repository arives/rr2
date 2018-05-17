#' Calculate R2.resid
#'
#' Calculate R2.resid for LMM, GLMM, PGLM, and PGLMMs.
#' @param mod a regression model with the following class: 'lmerMod', 'glmerMod', 'phylolm', and 'binaryPGLMM'
#' @param mod.r reduced model, if not provided, will use corresponding models with intercept as the only predictor
#' @param phy a phylogeny with tip labels and branch length
#' @return R2.resid
#' @export
#'
R2.resid <- function(mod = NULL, mod.r = NULL, phy = NULL) {
  
  if (!is.element(class(mod)[1], c("lm", "glm", "lmerMod", "glmerMod", "phylolm", "binaryPGLMM"))) {
    stop("mod must be class one of classes lm, glm, lmerMod, glmerMod, phylolm, binaryPGLMM.")
  }
  
  if (class(mod)[1] == "lm") {
    if (!is.object(mod.r)) {
      Y <- model.frame(mod)[, 1]
      mod.r <- lm(Y ~ 1)
    }
    if (!is.element(class(mod.r)[1], c("lm"))) {
      stop("mod.r must be class lm.")
    }
    return(R2.resid.lm(mod, mod.r))
  }
  
  if (class(mod)[1] == "glm") {
    if (!is.object(mod.r)) {
      Y <- model.frame(mod)[, 1]
      mod.r <- glm(Y ~ 1, family = family(mod)[[1]])
    }
    if (!is.element(class(mod.r)[1], c("glm"))) {
      stop("mod.r must be class glm.")
    }
    if (family(mod)[[1]] != family(mod.r)[[1]]) {
      stop("Sorry, but mod and mod.r must be from the same family of distributions.")
    }
    return(R2.resid.glm(mod, mod.r))
  }
  
  if (class(mod)[1] == "lmerMod") {
    if (!is.object(mod.r)) {
      # exists()?
      Y <- model.frame(mod)[, 1]
      mod.r <- lm(Y ~ 1)
    }
    if (!is.element(class(mod)[1], c("lmerMod", "lm"))) {
      stop("mod.r must be class lmerMod or lm.")
    }
    return(R2.resid.lmerMod(mod, mod.r))
  }
  
  if (class(mod)[1] == "glmerMod") {
    if (!is.object(mod.r)) {
      Y <- model.frame(mod)[, 1]
      mod.r <- glm(Y ~ 1, family = family(mod)[[1]])
    }
    if (!is.element(class(mod.r)[1], c("glmerMod", "glm"))) {
      stop("mod.r must be class glmerMod or glm.")
    }
    if (family(mod)[[1]] != family(mod.r)[[1]]) {
      stop("Sorry, but mod and mod.r must be from the same family of distributions.")
    }
    return(R2.resid.glmerMod(mod, mod.r))
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
    return(R2.resid.phylolm(mod, mod.r, phy))
  }
  
  if (class(mod)[1] == "binaryPGLMM") {
    if (!is.object(mod.r)) {
      Y <- mod$y
      mod.r <- glm(Y ~ 1, family = "binomial")
    }
    if (!is.element(class(mod.r)[1], c("binaryPGLMM", "glm"))) {
      stop("mod.r must be class binaryPGLMM or glm.")
    }
    return(R2.resid.binaryPGLMM(mod, mod.r))
  }
}

R2.resid.lm <- function(mod = NA, mod.r = NA) {
  X <- model.matrix(mod)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  X.r <- model.matrix(mod.r)
  p.r <- dim(X.r)[2]
  
  sigma2 <- (n - p)/n * stats::sigma(mod)^2
  sigma2.r <- (n - p.r)/n * stats::sigma(mod.r)^2
  R2.resid <- 1 - sigma2/sigma2.r
  return(R2.resid)
}

R2.resid.glm <- function(mod = NA, mod.r = NA) {
  mu <- mod$fitted.values
  if (family(mod)[1] == "binomial") 
    Yhat <- log(mu/(1 - mu))
  if (family(mod)[1] == "poisson") 
    Yhat <- log(mu)
  
  sig2e <- sig2e.r <- pi^2/3
  SSE.resid <- sig2e/(var(Yhat) + sig2e)
  
  mu.r <- mod.r$fitted.values
  if (family(mod.r)[1] == "binomial") 
    Yhat.r <- log(mu.r/(1 - mu.r))
  if (family(mod.r)[1] == "poisson") 
    Yhat.r <- log(mu.r)
  
  SSE.resid.r <- sig2e.r/(var(Yhat.r) + sig2e.r)
  
  R2.resid <- 1 - SSE.resid/SSE.resid.r
  return(R2.resid[1])
}

R2.resid.lmerMod <- function(mod = NULL, mod.r = NULL) {
  X <- model.matrix(mod)
  n <- dim(X)[1]
  X.r <- model.matrix(mod.r)
  p.r <- dim(X.r)[2]
  
  vcov <- as.data.frame(lme4::VarCorr(mod))$vcov
  sigma2 <- vcov[length(vcov)]
  
  if (class(mod.r) == "lmerMod") {
    vcov.r <- as.data.frame(lme4::VarCorr(mod.r))$vcov
    sigma2.r <- vcov.r[length(vcov.r)]
  }
  
  if (class(mod.r) == "lm") {
    sigma2.r <- (n - p.r)/n * stats::sigma(mod.r)^2
  }
  
  R2.resid <- 1 - sigma2/sigma2.r
  return(R2.resid)
}

R2.resid.glmerMod <- function(mod = NULL, mod.r = NULL) {
  
  X <- model.matrix(mod)
  n <- dim(X)[1]
  X.r <- model.matrix(mod.r)
  
  # full model
  AVAg <- crossprod(attr(mod, "pp")$LamtUt)
  Ag <- diag(1/attr(mod, "pp")$Xwts)
  C <- Ag %*% AVAg %*% Ag
  
  mu <- family(mod)$linkinv(X %*% lme4::fixef(mod))
  
  if (family(mod)[1] == "binomial") 
    Yhat <- log(mu/(1 - mu))
  if (family(mod)[1] == "poisson") 
    Yhat <- log(mu)
  
  sig2e <- sig2e.r <- pi^2/3
  sig2a <- prod(diag(C))^(1/n)
  
  SSE.resid <- sig2e/(var(Yhat) + sig2a + sig2e)
  
  # reduced model
  if (class(mod.r)[1] == "glmerMod") {
    AVAg.r <- crossprod(attr(mod.r, "pp")$LamtUt)
    Ag.r <- diag(1/attr(mod.r, "pp")$Xwts)
    
    C.r <- Ag.r %*% AVAg.r %*% Ag.r
    
    mu.r <- family(mod.r)$linkinv(X.r %*% lme4::fixef(mod.r))
    
    if (family(mod.r)[1] == "binomial") 
      Yhat.r <- log(mu.r/(1 - mu.r))
    if (family(mod.r)[1] == "poisson") 
      Yhat.r <- log(mu.r)
    
    sig2a.r <- prod(diag(C.r))^(1/n)
    
    SSE.resid.r <- sig2e.r/(var(Yhat.r) + sig2a.r + sig2e.r)
  }
  
  if (class(mod.r)[1] == "glm") {
    mu.r <- mod.r$fitted.values
    
    if (family(mod.r)[1] == "binomial") 
      Yhat.r <- log(mu.r/(1 - mu.r))
    if (family(mod.r)[1] == "poisson") 
      Yhat.r <- log(mu.r)
    
    SSE.resid.r <- sig2e.r/(var(Yhat.r) + sig2e.r)
  }
  
  R2.resid <- 1 - SSE.resid/SSE.resid.r
  return(R2.resid[1])
}

R2.resid.phylolm <- function(mod = NULL, mod.r = NULL, phy = NULL) {
  
  X <- mod$X
  n <- dim(X)[1]
  
  if (!mod$model %in% c("lambda", "OUrandomRoot", "OUfixedRoot", "BM", "kappa", "delta", "EB", "trend")) {
    stop("evolution model not supported yet")
  }
  
  phy.f <- transf_phy(mod, phy)
  
  scal <- sum(phy.f$edge.length)/n
  sigma2 <- mod$sigma2
  
  if (class(mod.r) == "phylolm") {
    if (!mod.r$model %in% c("lambda", "OUrandomRoot", "OUfixedRoot", "BM", "kappa", "delta", "EB", "trend")) {
      stop("evolution model not supported yet")
    }
    
    X.r <- mod.r$X
    p.r <- dim(X.r)[2]
    
    phy.r <- transf_phy(mod.r, phy)
    
    scal.r <- sum(phy.r$edge.length)/n
    sigma2.r <- mod.r$sigma2
  }
  
  if (class(mod.r) == "lm") {
    X.r <- model.matrix(mod.r)
    p.r <- dim(X.r)[2]
    scal.r <- 1
    sigma2.r <- (n - p.r)/n * stats::sigma(mod.r)^2
  }
  
  R2.resid <- 1 - (scal * sigma2)/(scal.r * sigma2.r)
  return(R2.resid)
}

R2.resid.binaryPGLMM <- function(mod = NULL, mod.r = NULL) {
  
  y <- mod$y
  n <- length(y)
  Yhat <- mod$X %*% mod$B
  phyV <- mod$VCV
  s2 <- mod$s2
  scal <- prod(diag(s2 * phyV))^(1/n)
  
  sig2e <- sig2e.r <- pi^2/3
  
  SSE.resid <- sig2e/(var(Yhat) + scal + sig2e)
  
  # reduced model
  if (class(mod.r)[1] == "binaryPGLMM") {
    Yhat.r <- mod.r$X %*% mod.r$B
    phyV.r <- mod.r$VCV
    s2.r <- mod.r$s2
    scal.r <- prod(diag(s2.r*phyV.r))^(1/n)
		
    SSE.resid.r <- sig2e.r/(var(Yhat.r) + scal.r + sig2e.r)
  }
  
  if (class(mod.r)[1] == "glm") {
    mu.r <- mod.r$fitted.values
    Yhat.r <- log(mu.r/(1 - mu.r))
    SSE.resid.r <- sig2e.r/(var(Yhat.r) + sig2e.r)
  }
  
  R2.resid <- 1 - SSE.resid/SSE.resid.r
  
  return(R2.resid[1])
}
