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
  
  if (!is.element(class(mod)[1], c("lm", "glm", "lmerMod", "glmerMod", "phylolm", "binaryPGLMM"))) {
    stop("mod must be class one of classes lm, glm, lmerMod, glmerMod, phylolm (but not phyloglm), binaryPGLMM.")
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
  Y <- model.frame(mod)[, 1]
  SSE.pred <- var(Y - lme4::predict(mod, re.form=NA))
  SSE.pred.r <- var(Y - lme4::predict(mod.r, re.form=NA))
  R2.pred <- 1 - SSE.pred/SSE.pred.r
  return(R2.pred)
}

R2.pred.glmerMod <- R2.pred.lmerMod

R2.pred.phylolm <- function(mod = NULL, mod.r = NULL, phy = NULL) {
  Y <- mod$y
  X <- mod$X
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  if (!mod$model %in% c("lambda", "OUrandomRoot", "OUfixedRoot", "BM", "kappa", "delta", "EB", "trend")) {
    stop("evolution model not supported yet")
  }
  
  if (!mod$model %in% c("BM", "trend")) {
    # optpar for BM models is NULL
    optpar <- round(mod$optpar, digits = 4)
    m.list <- list(x = optpar)
    if (mod$model %in% c("OUrandomRoot", "OUfixedRoot")) {
      names(m.list) <- "alpha"
    } else {
      names(m.list) <- mod$model
    }
    phy.f <- phylolm::transf.branch.lengths(phy, parameters = m.list, model = mod$model)$tree
  } else {
    phy.f <- phylolm::transf.branch.lengths(phy, parameters = NULL, model = mod$model)$tree
  }
  
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
    X.r <- mod.r$X
    if (!mod.r$model %in% c("lambda", "OUrandomRoot", "OUfixedRoot", "BM", "kappa", "delta", "EB", "trend")) {
      stop("evolution model not supported yet")
    }
    
    if (!mod.r$model %in% c("BM", "trend")) {
      optpar.r <- round(mod$optpar.r, digits = 4)
      m.list.r <- list(x = optpar.r)
      if (mod.r$model %in% c("OUrandomRoot", "OUfixedRoot")) {
        names(m.list.r) <- "alpha"
      } else {
        names(m.list.r) <- mod.r$model
      }
      phy.r <- phylolm::transf.branch.lengths(phy, parameters = m.list.r, model = mod.r$model)$tree
    } else {
      phy.r <- phylolm::transf.branch.lengths(phy, parameters = NULL, model = mod.r$model)$tree
    }
    
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
