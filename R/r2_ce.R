#' Calculate R2.ce
#'
#' Calculate R2.ce for LMM, GLMM, PGLM, and PGLMMs.
#' @param mod a regression model with the following class: "lmerMod", "glmerMod", "phylolm", and "binaryPGLMM"
#' @param mod.r reduced model, if not provided, will use corresponding models with intercept as the only predictor
#' @param phy a phylogeny with tip labels and branch length
#' @param cpp whether to use Rcpp version of loops, default is TRUE
#' @return R2.ce
#' @export
#'
R2.ce <- function(mod = NULL, mod.r = NULL, phy = NULL, cpp = TRUE) {

  if (!is.element(class(mod)[1], c("lmerMod", "glmerMod", "phylolm", "binaryPGLMM"))) {
    stop("mod must be class one of classes lmerMod, glmerMod, phylolm, binaryPGLMM.")
  }

  if (class(mod)[1] == "lmerMod") {
    if (!exists(deparse(substitute(mod.r)))) {
      Y <- model.frame(mod)[, 1]
      mod.r <- lm(Y ~ 1)
    }
    if (!is.element(class(mod)[1], c("lmerMod", "lm"))) {
      stop("mod.r must be class lmerMod or lm.")
    }
    if(cpp){
      return(R2.ce.lmerMod.cpp(mod, mod.r))
    } else {
      return(R2.ce.lmerMod(mod, mod.r))
    }
  }

  if (class(mod)[1] == "glmerMod") {
    if (family(mod)[[1]] != "binomial") {
      stop("Sorry, but only binomial (binary) glmerMod models are allowed.")
    }
    if (!exists(deparse(substitute(mod.r)))) {
      Y <- model.frame(mod)[, 1]
      mod.r <- glm(Y ~ 1, family = family(mod)[[1]])
    }
    if (!is.element(class(mod.r)[1], c("glmerMod", "glm"))) {
      stop("mod.r must be class glmerMod or glm.")
    }
    if(cpp){
      return(R2.ce.glmerMod.cpp(mod, mod.r)[1])
    } else {
      return(R2.ce.glmerMod(mod, mod.r)[1])
    }
  }

  if (class(mod)[1] == "phylolm") {
    if (!is.object(phy)) {
      stop("For phylolm you must provide the phylo object")
    }
    if (!exists(deparse(substitute(mod.r)))) {
      Y <- mod$y
      mod.r <- lm(Y ~ 1)
    }
    if (!is.element(class(mod.r)[1], c("phylolm", "lm"))) {
      stop("mod.r must be class phylolm or lm.")
    }
    if(cpp){
      return(R2.ce.phylolm.cpp(mod, mod.r, phy))
    } else {
      return(R2.ce.phylolm(mod, mod.r, phy))
    }
  }

  if (class(mod)[1] == "binaryPGLMM") {
    if (!exists(deparse(substitute(mod.r)))) {
      Y <- mod$y
      mod.r <- glm(Y ~ 1, family = "binomial")
    }
    if (!is.element(class(mod.r)[1], c("binaryPGLMM", "glm"))) {
      stop("mod.r must be class binaryPGLMM or glm.")
    }
    if(cpp){
      return(R2.ce.binaryPGLMM.cpp(mod, mod.r)[1])
    } else {
      return(R2.ce.binaryPGLMM(mod, mod.r)[1])
    }
  }
}

R2.ce.lmerMod <- function(mod = NULL, mod.r = NULL) {

  Y <- model.frame(mod)[, 1]
  X <- model.matrix(mod)
  n <- dim(X)[1]
  p <- dim(X)[2]

  X.r <- model.matrix(mod.r)
  p.r <- dim(X.r)[2]

  D <- attr(mod, "pp")$LamtUt
  V <- crossprod(D) + diag(n)

  if (class(mod.r) == "lmerMod") {
    D.r <- attr(mod.r, "pp")$LamtUt
    V.r <- crossprod(D.r) + diag(n)
  }
  if (class(mod.r) == "lm") {
    V.r <- diag(n)
  }

  iV <- solve(V)

  bhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% Y)
  fitted.values <- X %*% bhat
  R <- Y - X %*% bhat

  bbhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% R)[1] # Tony put this in the loop....

  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    VV <- V[-j, -j]
    iVV <- solve(VV)
    # version using global mean: This works best
    v <- V[j, -j]
    Rhat[j,] <- as.numeric(bbhat + v %*% iVV %*% (r - bbhat))
  }

  Yhat <- as.numeric(fitted.values + Rhat)
  SSE.ce <- var(Y - Yhat)

  # reduced model
  iV.r <- solve(V.r)

  bhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV.r %*% Y)
  fitted.values <- X.r %*% bhat
  R <- Y - X.r %*% bhat

  bbhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV %*% R)[1]

  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    VV.r <- V.r[-j, -j]
    iVV.r <- solve(VV.r)
    # version using global mean: This works best
    v <- V.r[j, -j]
    Rhat[j] <- as.numeric(bbhat + v %*% iVV.r %*% (r - bbhat))
  }

  Yhat.r <- as.numeric(fitted.values + Rhat)
  SSE.ce.r <- var(Y - Yhat.r)

  R2.ce <- 1 - SSE.ce/SSE.ce.r

  return(R2.ce)
}

R2.ce.lmerMod.cpp <- function(mod = NULL, mod.r = NULL) {

  Y <- model.frame(mod)[, 1]
  X <- model.matrix(mod)
  n <- dim(X)[1]
  p <- dim(X)[2]

  X.r <- model.matrix(mod.r)
  p.r <- dim(X.r)[2]

  D <- attr(mod, "pp")$LamtUt
  V <- crossprod(D) + diag(n)

  if (class(mod.r) == "lmerMod") {
    D.r <- attr(mod.r, "pp")$LamtUt
    V.r <- crossprod(D.r) + diag(n)
  }
  if (class(mod.r) == "lm") {
    V.r <- diag(n)
  }

  iV <- solve(V)

  bhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% Y)
  fitted.values <- X %*% bhat
  R <- Y - X %*% bhat

  Rhat = loop_cpp(as.matrix(R), as.matrix(V), as.matrix(iV), X)

  # Rhat <- matrix(0, nrow = n, ncol = 1)
  # for (j in 1:n) {
  #   r <- R[-j]
  #   VV <- V[-j, -j]
  #   iVV <- solve(VV)
  #   # version using global mean: This works best
  #   bbhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% R)[1]
  #   v <- V[j, -j]
  #   Rhat[j,] <- as.numeric(bbhat + v %*% iVV %*% (r - bbhat))
  # }

  Yhat <- as.numeric(fitted.values + Rhat)
  SSE.ce <- var(Y - Yhat)

  # reduced model
  iV.r <- solve(V.r)

  bhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV.r %*% Y)
  fitted.values <- X.r %*% bhat
  R <- Y - X.r %*% bhat

  Rhat = loop_cpp(as.matrix(R), as.matrix(V.r), as.matrix(iV.r), X.r)

  # Rhat <- matrix(0, nrow = n, ncol = 1)
  # for (j in 1:n) {
  #   r <- R[-j]
  #   VV.r <- V.r[-j, -j]
  #   iVV.r <- solve(VV.r)
  #   # version using global mean: This works best
  #   bbhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV %*% R)[1]
  #   v <- V.r[j, -j]
  #   Rhat[j] <- as.numeric(bbhat + v %*% iVV.r %*% (r - bbhat))
  # }

  Yhat.r <- as.numeric(fitted.values + Rhat)
  SSE.ce.r <- var(Y - Yhat.r)

  R2.ce <- 1 - SSE.ce/SSE.ce.r

  return(R2.ce)
}

R2.ce.glmerMod <- function(mod = NULL, mod.r = NULL) {

  Y <- model.frame(mod)[, 1]
  X <- model.matrix(mod)

  n <- dim(X)[1]
  p <- dim(X)[2]

  X.r <- model.matrix(mod.r)
  p.r <- dim(X.r)[2]

  AVAg <- crossprod(attr(mod, "pp")$LamtUt)
  Ag <- diag(1/attr(mod, "pp")$Xwts)

  mu <- family(mod)$linkinv(X %*% lme4::fixef(mod))
  mu[mu < 10^-6] <- 10^-6
  mu[mu > (1 - 10^-6)] <- (1 - 10^-6)

  g <- mu * (1 - mu)
  # g <- g/prod(g)^(1/n)
  V <- Ag %*% AVAg %*% Ag + diag(n)
  iV <- solve(V)

  iA <- diag(as.numeric(g^(-1/2)))
  iAVA <- iA %*% iV %*% iA

  A <- diag(as.numeric(g^(1/2)))
  AVA <- A %*% V %*% A

  Yhat <- mu

  # reduced model
  if (class(mod.r)[1] == "glmerMod") {
    AVAg.r <- crossprod(attr(mod.r, "pp")$LamtUt)
    Ag.r <- diag(1/attr(mod.r, "pp")$Xwts)

    mu.r <- family(mod.r)$linkinv(X.r %*% lme4::fixef(mod.r))
    mu.r[mu.r < 10^-6] <- 10^-6
    mu.r[mu.r > (1 - 10^-6)] <- (1 - 10^-6)

    g.r <- mu.r * (1 - mu.r)
    # g.r <- g.r/prod(g.r)^(1/n)
    V.r <- Ag.r %*% AVAg.r %*% Ag.r + diag(n)
    iV.r <- solve(V.r)

    iA.r <- diag(as.numeric(g.r^(-1/2)))
    iAVA.r <- iA.r %*% iV.r %*% iA.r

    A.r <- diag(as.numeric(g.r^(1/2)))
    AVA.r <- A.r %*% V.r %*% A.r

    Yhat.r <- mu.r
  }

  if (class(mod.r)[1] == "glm") {
    mu.r <- mod.r$fitted.values
    mu.r[mu.r < 10^-6] <- 10^-6
    mu.r[mu.r > (1 - 10^-6)] <- (1 - 10^-6)

    g.r <- mu.r * (1 - mu.r)
    # g.r <- g.r/prod(g.r)^(1/n)
    AVA.r <- diag(g.r)
    iAVA.r <- diag(1/g.r)
    V.r <- diag(n)
    iV.r <- diag(n)

    Yhat.r <- mu.r
  }

  R <- (Y - Yhat)
  X <- matrix(1, ncol = 1, nrow = n)
  # bbhat <- solve(t(X) %*% iAVA %*% X, t(X) %*% iAVA %*% R)[1]
  bbhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% R)[1]

  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    # VV <- AVA[-j,-j]
    VV <- V[-j, -j]
    iVV <- solve(VV)

    # v <- AVA[j, -j]
    v <- V[j, -j]
    Rhat[j] <- as.numeric(bbhat + v %*% iVV %*% (r - bbhat))
  }

  Yhat.f <- Yhat + Rhat
  SSE.ce <- var(Y - Yhat.f)

  # reduced model
  R.r <- (Y - Yhat.r)
  X <- matrix(1, ncol = 1, nrow = n)
  # bbhat.r <- solve(t(X) %*% iAVA.r %*% X, t(X) %*% iAVA.r %*% R.r)[1]
  bbhat.r <- solve(t(X) %*% iV.r %*% X, t(X) %*% iV.r %*% R.r)[1]

  Rhat.r <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r.r <- R.r[-j]
    # VV.r <- AVA.r[-j,-j]
    VV.r <- V.r[-j, -j]
    iVV.r <- solve(VV.r)

    # v.r <- AVA.r[j, -j]
    v.r <- V.r[j, -j]
    Rhat.r[j] <- as.numeric(bbhat.r + v.r %*% iVV.r %*% (r.r - bbhat.r))
  }

  Yhat.r <- Yhat.r + Rhat.r
  SSE.ce.r <- var(Y - Yhat.r)

  R2.ce <- 1 - SSE.ce/SSE.ce.r

  return(R2.ce)
}

R2.ce.glmerMod.cpp <- function(mod = NULL, mod.r = NULL) {

  Y <- model.frame(mod)[, 1]
  X <- model.matrix(mod)

  n <- dim(X)[1]
  p <- dim(X)[2]

  X.r <- model.matrix(mod.r)
  p.r <- dim(X.r)[2]

  AVAg <- crossprod(attr(mod, "pp")$LamtUt)
  Ag <- diag(1/attr(mod, "pp")$Xwts)

  mu <- family(mod)$linkinv(X %*% lme4::fixef(mod))
  mu[mu < 10^-6] <- 10^-6
  mu[mu > (1 - 10^-6)] <- (1 - 10^-6)

  g <- mu * (1 - mu)
  # g <- g/prod(g)^(1/n)
  V <- Ag %*% AVAg %*% Ag + diag(n)
  iV <- solve(V)

  iA <- diag(as.numeric(g^(-1/2)))
  iAVA <- iA %*% iV %*% iA

  A <- diag(as.numeric(g^(1/2)))
  AVA <- A %*% V %*% A

  Yhat <- mu

  # reduced model
  if (class(mod.r)[1] == "glmerMod") {
    AVAg.r <- crossprod(attr(mod.r, "pp")$LamtUt)
    Ag.r <- diag(1/attr(mod.r, "pp")$Xwts)

    mu.r <- family(mod.r)$linkinv(X.r %*% lme4::fixef(mod.r))
    mu.r[mu.r < 10^-6] <- 10^-6
    mu.r[mu.r > (1 - 10^-6)] <- (1 - 10^-6)

    g.r <- mu.r * (1 - mu.r)
    # g.r <- g.r/prod(g.r)^(1/n)
    V.r <- Ag.r %*% AVAg.r %*% Ag.r + diag(n)
    iV.r <- solve(V.r)

    iA.r <- diag(as.numeric(g.r^(-1/2)))
    iAVA.r <- iA.r %*% iV.r %*% iA.r

    A.r <- diag(as.numeric(g.r^(1/2)))
    AVA.r <- A.r %*% V.r %*% A.r

    Yhat.r <- mu.r
  }

  if (class(mod.r)[1] == "glm") {
    mu.r <- mod.r$fitted.values
    mu.r[mu.r < 10^-6] <- 10^-6
    mu.r[mu.r > (1 - 10^-6)] <- (1 - 10^-6)

    g.r <- mu.r * (1 - mu.r)
    # g.r <- g.r/prod(g.r)^(1/n)
    AVA.r <- diag(g.r)
    iAVA.r <- diag(1/g.r)
    V.r <- diag(n)
    iV.r <- diag(n)

    Yhat.r <- mu.r
  }

  R <- (Y - Yhat)
  X <- matrix(1, ncol = 1, nrow = n)

  Rhat = loop_cpp(as.matrix(R), as.matrix(V), as.matrix(iV), X)

  # # bbhat <- solve(t(X) %*% iAVA %*% X, t(X) %*% iAVA %*% R)[1]
  # bbhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% R)[1]
  # Rhat <- matrix(0, nrow = n, ncol = 1)
  # for (j in 1:n) {
  #   r <- R[-j]
  #   # VV <- AVA[-j,-j]
  #   VV <- V[-j, -j]
  #   iVV <- solve(VV)
  #
  #   # v <- AVA[j, -j]
  #   v <- V[j, -j]
  #   Rhat[j] <- as.numeric(bbhat + v %*% iVV %*% (r - bbhat))
  # }

  Yhat.f <- Yhat + Rhat
  SSE.ce <- var(Y - Yhat.f)

  # reduced model
  R.r <- (Y - Yhat.r)
  X <- matrix(1, ncol = 1, nrow = n)

  Rhat.r = loop_cpp(as.matrix(R.r), as.matrix(V.r), as.matrix(iV.r), X)

  # # bbhat.r <- solve(t(X) %*% iAVA.r %*% X, t(X) %*% iAVA.r %*% R.r)[1]
  # bbhat.r <- solve(t(X) %*% iV.r %*% X, t(X) %*% iV.r %*% R.r)[1]
  #
  # Rhat.r <- matrix(0, nrow = n, ncol = 1)
  # for (j in 1:n) {
  #   r.r <- R.r[-j]
  #   # VV.r <- AVA.r[-j,-j]
  #   VV.r <- V.r[-j, -j]
  #   iVV.r <- solve(VV.r)
  #
  #   # v.r <- AVA.r[j, -j]
  #   v.r <- V.r[j, -j]
  #   Rhat.r[j] <- as.numeric(bbhat.r + v.r %*% iVV.r %*% (r.r - bbhat.r))
  # }

  Yhat.r <- Yhat.r + Rhat.r
  SSE.ce.r <- var(Y - Yhat.r)

  R2.ce <- 1 - SSE.ce/SSE.ce.r

  return(R2.ce)
}

R2.ce.phylolm <- function(mod = NULL, mod.r = NULL, phy = NULL) {

  Y <- mod$y
  X <- mod$X
  n <- dim(X)[1]
  p <- dim(X)[2]

  optpar <- round(mod$optpar, digits = 4)
  if (mod$model == "lambda") {
    phy.f <- phylolm::transf.branch.lengths(phy,
                                            parameters = list(lambda = optpar),
                                            model = mod$model)$tree
  }
  if (mod$model %in% c("OUrandomRoot", "OUfixedRoot")){
    phy.f <- phylolm::transf.branch.lengths(phy,
                                            parameters = list(alpha = optpar),
                                            model = mod$model)$tree
  }

  V <- ape::vcv(phy.f)
  scal <- sum(phy.f$edge.length)/n
  V <- V/scal

  if (class(mod.r) == "phylolm") {
    X.r <- mod.r$X
    p.r <- dim(X.r)[2]

    optpar.r <- round(mod.r$optpar, digits = 4)

    if (mod$model == "lambda"){
      phy.r <- phylolm::transf.branch.lengths(phy,
                                              parameters = list(lambda = optpar.r),
                                              model = mod.r$model)$tree
    }

    if (mod$model %in% c("OUrandomRoot", "OUfixedRoot")){
      phy.r <- phylolm::transf.branch.lengths(phy,
                                              parameters = list(alpha = optpar.r),
                                              model = mod.r$model)$tree
    }

    V.r <- ape::vcv(phy.r)
    scal.r <- sum(phy.r$edge.length)/n
    V.r <- V.r/scal.r
  }

  if (class(mod.r) == "lm") {
    X.r <- model.matrix(mod.r)
    p.r <- dim(X.r)[2]
    V.r <- diag(n)
    scal.r <- 1
  }

  # full model
  iV <- solve(V)

  bhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% Y)
  fitted.values <- X %*% bhat
  R <- Y - X %*% bhat

  # version using global mean: This works best
  bbhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% R)[1]

  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    VV <- V[-j, -j]
    iVV <- solve(VV)
    v <- V[j, -j]
    Rhat[j] <- bbhat + v %*% iVV %*% (r - bbhat)
  }

  Yhat <- as.numeric(fitted.values + Rhat)
  SSE.ce <- var(Y - Yhat)

  # reduced model
  iV.r <- solve(V.r)

  bhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV.r %*% Y)
  fitted.values <- X.r %*% bhat
  R <- Y - X.r %*% bhat
  # version using global mean: This works best
  bbhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV %*% R)[1]

  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    VV.r <- V.r[-j, -j]
    iVV.r <- solve(VV.r)
    v <- V.r[j, -j]
    Rhat[j] <- bbhat + v %*% iVV.r %*% (r - bbhat)
  }

  Yhat.r <- as.numeric(fitted.values + Rhat)
  SSE.ce.r <- var(Y - Yhat.r)

  R2.ce <- 1 - SSE.ce/SSE.ce.r

  return(R2.ce)
}

R2.ce.phylolm.cpp <- function(mod = NULL, mod.r = NULL, phy = NULL) {

  Y <- mod$y
  X <- mod$X
  n <- dim(X)[1]
  p <- dim(X)[2]

  optpar <- round(mod$optpar, digits = 4)
  if (mod$model == "lambda") {
    phy.f <- phylolm::transf.branch.lengths(phy,
                                            parameters = list(lambda = optpar),
                                            model = mod$model)$tree
  }
  if (mod$model %in% c("OUrandomRoot", "OUfixedRoot")){
    phy.f <- phylolm::transf.branch.lengths(phy,
                                            parameters = list(alpha = optpar),
                                            model = mod$model)$tree
  }

  V <- ape::vcv(phy.f)
  scal <- sum(phy.f$edge.length)/n
  V <- V/scal

  if (class(mod.r) == "phylolm") {
    X.r <- mod.r$X
    p.r <- dim(X.r)[2]

    optpar.r <- round(mod.r$optpar, digits = 4)

    if (mod$model == "lambda"){
      phy.r <- phylolm::transf.branch.lengths(phy,
                                              parameters = list(lambda = optpar.r),
                                              model = mod.r$model)$tree
    }

    if (mod$model %in% c("OUrandomRoot", "OUfixedRoot")){
      phy.r <- phylolm::transf.branch.lengths(phy,
                                              parameters = list(alpha = optpar.r),
                                              model = mod.r$model)$tree
    }

    V.r <- ape::vcv(phy.r)
    scal.r <- sum(phy.r$edge.length)/n
    V.r <- V.r/scal.r
  }

  if (class(mod.r) == "lm") {
    X.r <- model.matrix(mod.r)
    p.r <- dim(X.r)[2]
    V.r <- diag(n)
    scal.r <- 1
  }

  # full model
  iV <- solve(V)

  bhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% Y)
  fitted.values <- X %*% bhat
  R <- Y - X %*% bhat

  Rhat = loop_cpp(as.matrix(R), as.matrix(V), as.matrix(iV), X)

  # Rhat <- matrix(0, nrow = n, ncol = 1)
  # for (j in 1:n) {
  #   r <- R[-j]
  #   VV <- V[-j, -j]
  #   iVV <- solve(VV)
  #
  #   # version using global mean: This works best
  #   bbhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% R)[1]
  #   v <- V[j, -j]
  #   Rhat[j] <- bbhat + v %*% iVV %*% (r - bbhat)
  # }

  Yhat <- as.numeric(fitted.values + Rhat)
  SSE.ce <- var(Y - Yhat)

  # reduced model
  iV.r <- solve(V.r)

  bhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV.r %*% Y)
  fitted.values <- X.r %*% bhat
  R <- Y - X.r %*% bhat

  Rhat = loop_cpp(as.matrix(R), as.matrix(V.r), as.matrix(iV.r), X.r)

  # Rhat <- matrix(0, nrow = n, ncol = 1)
  # for (j in 1:n) {
  #   r <- R[-j]
  #   VV.r <- V.r[-j, -j]
  #   iVV.r <- solve(VV.r)
  #   # version using global mean: This works best
  #   bbhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV %*% R)[1]
  #   v <- V.r[j, -j]
  #   Rhat[j] <- bbhat + v %*% iVV.r %*% (r - bbhat)
  # }

  Yhat.r <- as.numeric(fitted.values + Rhat)
  SSE.ce.r <- var(Y - Yhat.r)

  R2.ce <- 1 - SSE.ce/SSE.ce.r

  return(R2.ce)
}

R2.ce.binaryPGLMM <- function(mod = NULL, mod.r = NULL) {

  mu <- mod$mu
  n <- length(mu)
  Yhat <- mu
  Y <- mod$y
  R <- Y - mu

  if (class(mod.r)[1] == "binaryPGLMM") {
    mu.r <- mod.r$mu
    Yhat.r <- mu.r
    Y <- mod.r$y
    R.r <- Y - mu.r
  } else {
    mu.r <- mod.r$fitted.values
    Yhat.r <- mu.r
    Y <- mod.r$y
    R.r <- Y - mu.r
  }

  phy <- mod$phy

  g <- mu * (1 - mu)
  g <- g/prod(g)^(1/n)
  V <- mod$s2 * ape::vcv(phy) + diag(n)
  phy.temp <- ape::vcv2phylo(V)
  V <- V/(sum(phy.temp$edge.length)/n)
  iV <- solve(V)

  iA <- diag(as.numeric(g^(-1/2)))
  iAVA <- iA %*% iV %*% iA

  A <- diag(as.numeric(g^(1/2)))
  AVA <- A %*% V %*% A

  Yhat <- mu

  # reduced model
  if (class(mod.r)[1] == "binaryPGLMM") {
    g.r <- mu.r * (1 - mu.r)
    g.r <- g.r/prod(g.r)^(1/n)
    V.r <- mod.r$s2 * ape::vcv(phy) + diag(n)
    phy.temp <- ape::vcv2phylo(V.r)
    V.r <- V.r/(sum(phy.temp$edge.length)/n)
    iV.r <- solve(V.r)

    iA.r <- diag(as.numeric(g.r^(-1/2)))
    iAVA.r <- iA.r %*% iV.r %*% iA.r

    A.r <- diag(as.numeric(g.r^(1/2)))
    AVA.r <- A.r %*% V.r %*% A.r
  }

  if (class(mod.r)[1] == "glm") {
    mu.r <- mod.r$fitted.values
    g.r <- mu.r * (1 - mu.r)
    g.r <- g.r/prod(g.r)^(1/n)
    AVA.r <- diag(g.r)
    iAVA.r <- diag(1/g.r)
  }

  R <- (Y - Yhat)
  X <- matrix(1, ncol = 1, nrow = n)
  bbhat <- solve(t(X) %*% iAVA %*% X, t(X) %*% iAVA %*% R)[1]

  Rhat <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r <- R[-j]
    VV <- AVA[-j, -j]
    iVV <- solve(VV)

    v <- AVA[j, -j]
    Rhat[j] <- as.numeric(bbhat + v %*% iVV %*% (r - bbhat))
  }

  Yhat.f <- Yhat + Rhat
  SSE.ce <- var(Y - Yhat.f)

  # reduced model
  R.r <- (Y - Yhat.r)
  X <- matrix(1, ncol = 1, nrow = n)
  bbhat.r <- solve(t(X) %*% iAVA.r %*% X, t(X) %*% iAVA.r %*% R.r)[1]

  Rhat.r <- matrix(0, nrow = n, ncol = 1)
  for (j in 1:n) {
    r.r <- R.r[-j]
    VV.r <- AVA.r[-j, -j]
    iVV.r <- solve(VV.r)

    v.r <- AVA.r[j, -j]
    Rhat.r[j] <- as.numeric(bbhat.r + v.r %*% iVV.r %*% (r.r - bbhat.r))
  }

  Yhat.r <- Yhat.r + Rhat.r
  SSE.ce.r <- var(Y - Yhat.r)

  R2.ce <- 1 - SSE.ce/SSE.ce.r

  return(R2.ce)
}

R2.ce.binaryPGLMM.cpp <- function(mod = NULL, mod.r = NULL) {

  mu <- mod$mu
  n <- length(mu)
  Yhat <- mu
  Y <- mod$y
  R <- Y - mu

  if (class(mod.r)[1] == "binaryPGLMM") {
    mu.r <- mod.r$mu
    Yhat.r <- mu.r
    Y <- mod.r$y
    R.r <- Y - mu.r
  } else {
    mu.r <- mod.r$fitted.values
    Yhat.r <- mu.r
    Y <- mod.r$y
    R.r <- Y - mu.r
  }

  phy <- mod$phy

  g <- mu * (1 - mu)
  g <- g/prod(g)^(1/n)
  V <- mod$s2 * ape::vcv(phy) + diag(n)
  phy.temp <- ape::vcv2phylo(V)
  V <- V/(sum(phy.temp$edge.length)/n)
  iV <- solve(V)

  iA <- diag(as.numeric(g^(-1/2)))
  iAVA <- iA %*% iV %*% iA

  A <- diag(as.numeric(g^(1/2)))
  AVA <- A %*% V %*% A

  Yhat <- mu

  # reduced model
  if (class(mod.r)[1] == "binaryPGLMM") {
    g.r <- mu.r * (1 - mu.r)
    g.r <- g.r/prod(g.r)^(1/n)
    V.r <- mod.r$s2 * ape::vcv(phy) + diag(n)
    phy.temp <- ape::vcv2phylo(V.r)
    V.r <- V.r/(sum(phy.temp$edge.length)/n)
    iV.r <- solve(V.r)

    iA.r <- diag(as.numeric(g.r^(-1/2)))
    iAVA.r <- iA.r %*% iV.r %*% iA.r

    A.r <- diag(as.numeric(g.r^(1/2)))
    AVA.r <- A.r %*% V.r %*% A.r
  }

  if (class(mod.r)[1] == "glm") {
    mu.r <- mod.r$fitted.values
    g.r <- mu.r * (1 - mu.r)
    g.r <- g.r/prod(g.r)^(1/n)
    AVA.r <- diag(g.r)
    iAVA.r <- diag(1/g.r)
  }

  R <- (Y - Yhat)
  X <- matrix(1, ncol = 1, nrow = n)

  Rhat = loop_cpp(as.matrix(R), as.matrix(AVA), as.matrix(iAVA), X)
  # bbhat <- solve(t(X) %*% iAVA %*% X, t(X) %*% iAVA %*% R)[1]
  # Rhat <- matrix(0, nrow = n, ncol = 1)
  # for (j in 1:n) {
  #   r <- R[-j]
  #   VV <- AVA[-j, -j]
  #   iVV <- solve(VV)
  #
  #   v <- AVA[j, -j]
  #   Rhat[j] <- as.numeric(bbhat + v %*% iVV %*% (r - bbhat))
  # }

  Yhat.f <- Yhat + Rhat
  SSE.ce <- var(Y - Yhat.f)

  # reduced model
  R.r <- (Y - Yhat.r)
  X <- matrix(1, ncol = 1, nrow = n)

  Rhat.r = loop_cpp(as.matrix(R.r), as.matrix(AVA.r), as.matrix(iAVA.r), X)

  # bbhat.r <- solve(t(X) %*% iAVA.r %*% X, t(X) %*% iAVA.r %*% R.r)[1]
  # Rhat.r <- matrix(0, nrow = n, ncol = 1)
  # for (j in 1:n) {
  #   r.r <- R.r[-j]
  #   VV.r <- AVA.r[-j, -j]
  #   iVV.r <- solve(VV.r)
  #
  #   v.r <- AVA.r[j, -j]
  #   Rhat.r[j] <- as.numeric(bbhat.r + v.r %*% iVV.r %*% (r.r - bbhat.r))
  # }

  Yhat.r <- Yhat.r + Rhat.r
  SSE.ce.r <- var(Y - Yhat.r)

  R2.ce <- 1 - SSE.ce/SSE.ce.r

  return(R2.ce)
}
