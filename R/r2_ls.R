#' Calculate R2.ls
#'
#' Calculate R2.ls for LMM, GLMM, PGLM, and PGLMMs.
#' @param mod a regression model with the following class: "lmerMod", "glmerMod", "phylolm", and "binaryPGLMM"
#' @param mod.r reduced model, if not provided, will use corresponding models with intercept as the only predictor
#' @param phy a phylogeny with tip labels and branch length
#' @return R2.ls
#' @export
#'
R2.ls <- function(mod = NULL, mod.r = NULL, phy = NULL) {

  if (!is.element(class(mod)[1], c("lmerMod", "glmerMod", "phylolm", "binaryPGLMM"))) {
    stop("mod must be class one of classes lmerMod, glmerMod, phylolm, binaryPGLMM.")
  }

  if (class(mod)[1] == "lmerMod") {
    if (!exists(deparse(substitute(mod.r)))) {
      # exists()?
      Y <- model.frame(mod)[, 1]
      mod.r <- lm(Y ~ 1)
    }
    if (!is.element(class(mod)[1], c("lmerMod", "lm"))) {
      stop("mod.r must be class lmerMod or lm.")
    }
    return(R2.ls.lmerMod(mod, mod.r))
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
    return(R2.ls.glmerMod(mod, mod.r)[1])
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
    return(R2.ls.phylolm(mod, mod.r, phy))
  }

  if (class(mod)[1] == "binaryPGLMM") {
    if (!exists(deparse(substitute(mod.r)))) {
      Y <- mod$y
      mod.r <- glm(Y ~ 1, family = "binomial")
    }
    if (!is.element(class(mod.r)[1], c("binaryPGLMM", "glm"))) {
      stop("mod.r must be class binaryPGLMM or glm.")
    }
    return(R2.ls.binaryPGLMM(mod, mod.r)[1])
  }
}

R2.ls.lmerMod <- function(mod = NULL, mod.r = NULL) {

  Y <- model.frame(mod)[,1]
  X <- model.matrix(mod)
  n <- dim(X)[1]
  p <- dim(X)[2]

  X.r <- model.matrix(mod.r)
  p.r <- dim(X.r)[2]

  D <- attr(mod,"pp")$LamtUt
  V <- crossprod(D) + diag(n) # crossprod(D) is faster than t(D) %*% D

  if(class(mod.r) == "lmerMod"){
    D.r <- attr(mod.r,"pp")$LamtUt
    V.r <- crossprod(D.r) + diag(n)
  }

  if(class(mod.r) == "lm"){
    V.r <- diag(n)
  }

  vcov <- as.data.frame(lme4::VarCorr(mod))$vcov

  sigma2 <- vcov[length(vcov)]

  if(class(mod.r) == "lmerMod"){
    vcov.r <- as.data.frame(lme4::VarCorr(mod.r))$vcov
    sigma2.r <- vcov.r[length(vcov.r)]
  }

  if(class(mod.r) == "lm"){
    sigma2.r <- (n - p.r)/n * stats::sigma(mod.r)^2
  }

  R2.ls <- 1 - sigma2/sigma2.r

  return(R2.ls)
}

R2.ls.glmerMod <- function(mod = NULL, mod.r = NULL) {

  Y <- model.frame(mod)[, 1]
  X <- model.matrix(mod)
  n <- dim(X)[1]
  p <- dim(X)[2]

  X.r <- model.matrix(mod.r)
  p.r <- dim(X.r)[2]

  # full model
  AVAg <- crossprod(attr(mod, "pp")$LamtUt)
  Ag <- diag(1/attr(mod, "pp")$Xwts)

  mu <- family(mod)$linkinv(X %*% lme4::fixef(mod))

  g <- mu * (1 - mu)
  g <- g/prod(g)^(1/n)
  V <- Ag %*% AVAg %*% Ag + diag(n)
  iV <- solve(V)

  iA <- diag(as.numeric(g^(-1/2)))
  iAVA <- iA %*% iV %*% iA
  SSE.ls <- t(Y - mu) %*% iAVA %*% (Y - mu)

  # reduced model
  if (class(mod.r)[1] == "glmerMod") {
    AVAg.r <- crossprod(attr(mod.r, "pp")$LamtUt)
    Ag.r <- diag(1/attr(mod.r, "pp")$Xwts)

    mu.r <- family(mod.r)$linkinv(X.r %*% lme4::fixef(mod.r))
    g.r <- mu.r * (1 - mu.r)
    g.r <- g.r/prod(g.r)^(1/n)
    V.r <- Ag.r %*% AVAg.r %*% Ag.r + diag(n)
    iV.r <- solve(V.r)

    iA.r <- diag(as.numeric(g.r^(-1/2)))
    iAVA.r <- iA.r %*% iV.r %*% iA.r
    SSE.ls.r <- t(Y - mu.r) %*% iAVA.r %*% (Y - mu.r)
  }

  if (class(mod.r)[1] == "glm") {
    mu.r <- mod.r$fitted.values
    g.r <- mu.r * (1 - mu.r)
    g.r <- g.r/prod(g.r)^(1/n)
    V.r <- diag(g.r)
    iV.r <- diag(1/g.r)
    SSE.ls.r <- t(Y - mu.r) %*% iV.r %*% (Y - mu.r)
  }

  R2.ls <- 1 - SSE.ls/SSE.ls.r

  return(R2.ls)
}

R2.ls.phylolm <- function(mod = NULL, mod.r = NULL, phy = NULL) {

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

  if (mod$model %in% c("OUrandomRoot", "OUfixedRoot")) {
    phy.f <- phylolm::transf.branch.lengths(phy,
                                            parameters = list(alpha = optpar),
                                            model = mod$model)$tree
  }

  V <- ape::vcv(phy.f)
  scal <- sum(phy.f$edge.length)/n
  V <- V/scal

  sigma2 <- mod$sigma2

  if (class(mod.r) == "phylolm") {
    X.r <- mod.r$X
    p.r <- dim(X.r)[2]

    optpar.r <- round(mod.r$optpar, digits = 4)

    if (mod$model == "lambda") {
      phy.r <- phylolm::transf.branch.lengths(phy,
                                              parameters = list(lambda = optpar.r),
                                              model = mod.r$model)$tree
    }

    if (mod$model %in% c("OUrandomRoot", "OUfixedRoot")) {
      phy.r <- phylolm::transf.branch.lengths(phy,
                                              parameters = list(alpha = optpar.r),
                                              model = mod.r$model)$tree
    }

    V.r <- ape::vcv(phy.r)
    scal.r <- sum(phy.r$edge.length)/n
    V.r <- V.r/scal.r
    sigma2.r <- mod.r$sigma2
  }

  if (class(mod.r) == "lm") {
    X.r <- model.matrix(mod.r)
    p.r <- dim(X.r)[2]
    V.r <- diag(n)
    scal.r <- 1
    sigma2.r <- (n - p.r)/n * sigma(mod.r)^2
  }

  R2.ls <- 1 - (scal * sigma2)/(scal.r * sigma2.r)

  return(R2.ls)
}

R2.ls.binaryPGLMM <- function(mod = NULL, mod.r = NULL) {
  mu <- mod$mu
  n <- length(mu)
  Yhat <- mu
  R <- mod$y - mu

  if (class(mod.r)[1] == "binaryPGLMM") {
    mu.r <- mod.r$mu
    Yhat.r <- mu.r
    R.r <- mod.r$y - mu.r
  } else {
    mu.r <- mod.r$fitted.values
    Yhat.r <- mu.r
    R.r <- mod.r$y - mu.r
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
  SSE.ls <- t(R) %*% iAVA %*% R

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
    SSE.ls.r <- t(R.r) %*% iAVA.r %*% R.r
  }

  if (class(mod.r)[1] == "glm") {
    mu.r <- mod.r$fitted.values
    g.r <- mu.r * (1 - mu.r)
    g.r <- g.r/prod(g.r)^(1/n)
    V.r <- diag(g.r)
    iV.r <- diag(1/g.r)
    SSE.ls.r <- t(R.r) %*% iV.r %*% R.r
  }

  R2.ls <- 1 - SSE.ls/SSE.ls.r

  return(R2.ls)
}
