context("testing whether rr2 pkg get same results as original code")

test_that("results from rr2 package should be the same as the original R code", {
  #### begin of original code --------- original R functions, renamed by adding .original suff R2_source_code
  library(lme4)
  library(phylolm)
  library(ape)
  
  ######################################################## R2.ls
  
  R2.ls.original <- function(mod = NA, mod.r = NA, phy = NA) {
    
    if (!is.element(class(mod)[1], c("lmerMod", "glmerMod", "phylolm", "binaryPGLMM"))) {
      stop("mod must be class one of classes lmerMod, glmerMod, phylolm, binaryPGLMM.")
    }
    if (class(mod)[1] == "lmerMod") {
      if (!is.object(mod.r)) {
        Y <- model.frame(mod)[, 1]
        mod.r <- lm(Y ~ 1)
      }
      if (!is.element(class(mod)[1], c("lmerMod", "lm"))) {
        stop("mod.r must be class lmerMod or lm.")
      }
      return(R2.ls.lmerMod.original(mod, mod.r))
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
      return(R2.ls.glmerMod.original(mod, mod.r))
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
      return(R2.ls.phylolm.original(mod, mod.r, phy))
    }
    if (class(mod)[1] == "binaryPGLMM") {
      if (!is.object(mod.r)) {
        Y <- mod$y
        mod.r <- glm(Y ~ 1, family = "binomial")
      }
      if (!is.element(class(mod.r)[1], c("binaryPGLMM", "glm"))) {
        stop("mod.r must be class binaryPGLMM or glm.")
      }
      return(R2.ls.binaryPGLMM.original(mod, mod.r))
    }
  }
  
  R2.ls.lmerMod.original <- function(mod = NA, mod.r = NA) {
    
    Y <- model.frame(mod)[, 1]
    X <- model.matrix(mod)
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    X.r <- model.matrix(mod.r)
    p.r <- dim(X.r)[2]
    
    vcov <- as.data.frame(VarCorr(mod))$vcov
    sigma2 <- vcov[length(vcov)]
    if (class(mod.r) == "lmerMod") {
      vcov.r <- as.data.frame(VarCorr(mod.r))$vcov
      sigma2.r <- vcov.r[length(vcov.r)]
    }
    if (class(mod.r) == "lm") {
      sigma2.r <- (n - p.r)/n * sigma(mod.r)^2
    }
    R2.ls <- 1 - sigma2/sigma2.r
    return(R2.ls)
  }
  
  R2.ls.glmerMod.original <- function(mod = NA, mod.r = NA) {
    
    Y <- model.frame(mod)[, 1]
    X <- model.matrix(mod)
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    X.r <- model.matrix(mod.r)
    p.r <- dim(X.r)[2]
    
    # full model
    AVAg <- t(attr(mod, "pp")$LamtUt) %*% attr(mod, "pp")$LamtUt
    Ag <- diag(1/attr(mod, "pp")$Xwts)
    C <- Ag %*% AVAg %*% Ag
    
    mu <- family(mod)$linkinv(X %*% fixef(mod))
    
    if (family(mod)[1] == "binomial") 
      v <- mu * (1 - mu)
    if (family(mod)[1] == "poisson") 
      v <- mu
    
    sig2e <- pi^2/3
    sig2a <- prod(diag(C))^(1/n)
    # sig2a <- mean(diag(C))
    Yhat <- log(mu/(1 - mu))
    
    SSE.ls <- sig2e/(var(Yhat) + sig2a + sig2e)
    
    # reduced model
    if (class(mod.r)[1] == "glmerMod") {
      AVAg.r <- t(attr(mod.r, "pp")$LamtUt) %*% attr(mod.r, "pp")$LamtUt
      Ag.r <- diag(1/attr(mod.r, "pp")$Xwts)
      
      C.r <- Ag.r %*% AVAg.r %*% Ag.r
      
      mu.r <- family(mod.r)$linkinv(X.r %*% fixef(mod.r))
      
      if (family(mod.r)[1] == "binomial") 
        v.r <- mu.r * (1 - mu.r)
      if (family(mod.r)[1] == "poisson") 
        v.r <- mu.r
      
      sig2e.r <- pi^2/3
      sig2a.r <- prod(diag(C.r))^(1/n)
      # sig2a.r <- mean(diag(C.r))
      Yhat.r <- log(mu.r/(1 - mu.r))
      
      SSE.ls.r <- sig2e.r/(var(Yhat.r) + sig2a.r + sig2e.r)
    }
    if (class(mod.r)[1] == "glm") {
      mu.r <- mod.r$fitted.values
      
      if (family(mod.r)[1] == "binomial") 
        v.r <- mu.r * (1 - mu.r)
      if (family(mod.r)[1] == "poisson") 
        v.r <- mu.r
      
      sig2e.r <- pi^2/3
      Yhat.r <- log(mu.r/(1 - mu.r))
      
      SSE.ls.r <- sig2e.r/(var(Yhat.r) + sig2e.r)
    }
    R2.ls <- 1 - SSE.ls/SSE.ls.r
    return(R2.ls[1])
  }
  
  R2.ls.phylolm.original <- function(mod = NA, mod.r = NA, phy = NA) {
    
    Y <- mod$y
    X <- mod$X
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    optpar <- round(mod$optpar, digits = 4)
    if (mod$model == "lambda") 
      phy.f <- transf.branch.lengths(phy, parameters = list(lambda = optpar), model = mod$model)$tree
    if (mod$model == "OUrandomRoot") 
      phy.f <- transf.branch.lengths(phy, parameters = list(alpha = optpar), model = mod$model)$tree
    if (mod$model == "OUfixedRoot") 
      phy.f <- transf.branch.lengths(phy, parameters = list(alpha = optpar), model = mod$model)$tree
    
    scal <- sum(phy.f$edge.length)/n
    sigma2 <- mod$sigma2
    
    if (class(mod.r) == "phylolm") {
      X.r <- mod.r$X
      p.r <- dim(X.r)[2]
      
      optpar.r <- round(mod.r$optpar, digits = 4)
      if (mod.r$model == "lambda") 
        phy.r <- transf.branch.lengths(phy, parameters = list(lambda = optpar.r), model = mod.r$model)$tree
      if (mod.r$model == "OUrandomRoot") 
        phy.r <- transf.branch.lengths(phy, parameters = list(alpha = optpar.r), model = mod.r$model)$tree
      if (mod.r$model == "OUfixedRoot") 
        phy.r <- transf.branch.lengths(phy, parameters = list(alpha = optpar.r), model = mod.r$model)$tree
      
      scal.r <- sum(phy.r$edge.length)/n
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
  
  R2.ls.binaryPGLMM.original <- function(mod = NA, mod.r = NA) {
    
    y <- mod$y
    n <- length(y)
    Yhat <- mod$X %*% mod$B
    phyV <- mod$VCV
    s2 <- mod$s2
    
    sig2e <- pi^2/3
    
    SSE.ls <- sig2e/(var(Yhat) + s2 + sig2e)
    
    # reduced model
    if (class(mod.r)[1] == "binaryPGLMM") {
      Yhat.r <- mod.r$X %*% mod.r$B
      phyV.r <- mod.r$VCV
      s2.r <- mod.r$s2
      
      sig2e.r <- pi^2/3
      
      SSE.ls.r <- sig2e.r/(var(Yhat.r) + s2.r + sig2e.r)
    }
    if (class(mod.r)[1] == "glm") {
      mu.r <- mod.r$fitted.values
      
      sig2e.r <- pi^2/3
      Yhat.r <- log(mu.r/(1 - mu.r))
      
      SSE.ls.r <- sig2e.r/(var(Yhat.r) + sig2e.r)
    }
    
    R2.ls <- 1 - SSE.ls/SSE.ls.r
    
    return(R2.ls[1])
  }
  
  
  ######################################################## R2.ce
  
  R2.ce.original <- function(mod = NA, mod.r = NA, phy = NA) {
    
    if (!is.element(class(mod)[1], c("lmerMod", "glmerMod", "phylolm", "binaryPGLMM"))) {
      stop("mod must be class one of classes lmerMod, glmerMod, phylolm (but not phyloglm), binaryPGLMM.")
    }
    if (class(mod)[1] == "lmerMod") {
      if (!is.object(mod.r)) {
        Y <- model.frame(mod)[, 1]
        mod.r <- lm(Y ~ 1)
      }
      if (!is.element(class(mod)[1], c("lmerMod", "lm"))) {
        stop("mod.r must be class lmerMod or lm.")
      }
      return(R2.ce.lmerMod.original(mod, mod.r))
    }
    if (class(mod)[1] == "glmerMod") {
      if (!is.object(mod.r)) {
        Y <- model.frame(mod)[, 1]
        mod.r <- glm(Y ~ 1, family = family(mod)[[1]])
      }
      if (!is.element(class(mod.r)[1], c("glmerMod", "glm"))) {
        stop("mod.r must be class glmerMod or glm.")
      }
      return(R2.ce.glmerMod.original(mod, mod.r))
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
      return(R2.ce.phylolm.original(mod, mod.r, phy))
    }
    if (class(mod)[1] == "binaryPGLMM") {
      if (!is.object(mod.r)) {
        Y <- mod$y
        mod.r <- glm(Y ~ 1, family = "binomial")
      }
      if (!is.element(class(mod.r)[1], c("binaryPGLMM", "glm"))) {
        stop("mod.r must be class binaryPGLMM or glm.")
      }
      return(R2.ce.binaryPGLMM.original(mod, mod.r))
    }
  }
  
  R2.ce.lmerMod.original <- function(mod = NA, mod.r = NA) {
    
    Y <- model.frame(mod)[, 1]
    SSE.ce <- var(Y - fitted(mod))
    
    # reduced model
    if (class(mod.r)[1] == "lmerMod") {
      SSE.ce.r <- var(Y - fitted(mod.r))
    }
    if (class(mod.r)[1] == "lm") {
      SSE.ce.r <- var(Y - mod.r$fitted.values)
    }
    
    R2.ce <- 1 - SSE.ce/SSE.ce.r
    return(R2.ce)
  }
  
  R2.ce.glmerMod.original <- function(mod = NA, mod.r = NA) {
    
    Y <- model.frame(mod)[, 1]
    SSE.ce <- var(Y - fitted(mod))
    SSE.ce.r <- var(Y - fitted(mod.r))
    R2.ce <- 1 - SSE.ce/SSE.ce.r
    return(R2.ce)
  }
  
  R2.ce.phylolm.original <- function(mod = NA, mod.r = NA, phy = NA) {
    
    Y <- mod$y
    X <- mod$X
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    optpar <- round(mod$optpar, digits = 4)
    if (mod$model == "lambda") 
      phy.f <- transf.branch.lengths(phy, parameters = list(lambda = optpar), model = mod$model)$tree
    if (mod$model == "OUrandomRoot") 
      phy.f <- transf.branch.lengths(phy, parameters = list(alpha = optpar), model = mod$model)$tree
    if (mod$model == "OUfixedRoot") 
      phy.f <- transf.branch.lengths(phy, parameters = list(alpha = optpar), model = mod$model)$tree
    V <- vcv(phy.f)
    
    R <- Y - fitted(mod)
    Rhat <- matrix(0, nrow = n, ncol = 1)
    for (j in 1:n) {
      
      r <- R[-j]
      VV <- V[-j, -j]
      iVV <- solve(VV)
      
      v <- V[j, -j]
      Rhat[j] <- v %*% iVV %*% r
    }
    Yhat <- as.numeric(fitted(mod) + Rhat)
    SSE.ce <- var(Y - Yhat)
    
    # reduced model
    if (class(mod.r) == "phylolm") {
      X.r <- mod.r$X
      
      optpar.r <- round(mod.r$optpar, digits = 4)
      if (mod$model == "lambda") 
        phy.r <- transf.branch.lengths(phy, parameters = list(lambda = optpar.r), model = mod.r$model)$tree
      if (mod$model == "OUrandomRoot") 
        phy.r <- transf.branch.lengths(phy, parameters = list(alpha = optpar.r), model = mod.r$model)$tree
      if (mod$model == "OUfixedRoot") 
        phy.r <- transf.branch.lengths(phy, parameters = list(alpha = optpar.r), model = mod.r$model)$tree
      V.r <- vcv(phy.r)
      
      R.r <- Y - fitted(mod.r)
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
      Yhat.r <- fitted(mod.r)
    }
    SSE.ce.r <- var(Y - Yhat.r)
    
    R2.ce <- 1 - SSE.ce/SSE.ce.r
    return(R2.ce[1])
  }
  
  R2.ce.binaryPGLMM.original <- function(mod = NA, mod.r = NA) {
    
    Yhat <- mod$mu
    Y <- mod$y
    
    SSE.ce <- var(Y - Yhat)
    
    if (class(mod.r)[1] == "binaryPGLMM") {
      Yhat.r <- mod.r$mu
      SSE.ce.r <- var(Y - Yhat.r)
    }
    if (class(mod.r)[1] == "glm") {
      Yhat.r <- mod.r$fitted.values
      SSE.ce.r <- var(Y - Yhat.r)
    }
    
    R2.ce <- 1 - SSE.ce/SSE.ce.r
    return(R2.ce[1])
  }
  
  ######################################################## R2.lr
  
  R2.lr.original <- function(mod = NA, mod.r = NA) {
    
    if (!is.element(class(mod)[1], c("lmerMod", "glmerMod", "phylolm", "phyloglm"))) {
      stop("mod must be class one of classes lmerMod, glmerMod, phylolm, phyloglm.")
    }
    if (class(mod)[1] == "lmerMod") {
      if (!is.object(mod.r)) {
        Y <- model.frame(mod)[, 1]
        mod.r <- lm(Y ~ 1)
      }
      if (!is.element(class(mod)[1], c("lmerMod", "lm"))) {
        stop("mod.r must be class lmerMod or lm.")
      }
      if (isREML(mod)) {
        mod <- update(mod, REML = F)
        warning("mod updated with REML=F")
      }
      if (class(mod.r)[1] == "lmerMod" && isREML(mod.r)) {
        mod.r <- update(mod.r, REML = F)
        warning("mod.r updated with REML=F")
      }
      return(R2.lr.lmerMod.original(mod, mod.r)[1])
    }
    if (class(mod)[1] == "glmerMod") {
      if (!is.object(mod.r)) {
        Y <- model.frame(mod)[, 1]
        mod.r <- glm(Y ~ 1, family = family(mod)[[1]])
      }
      if (!is.element(class(mod.r)[1], c("glmerMod", "glm"))) {
        stop("mod.r must be class glmerMod or glm.")
      }
      return(R2.lr.glmerMod.original(mod, mod.r)[1])
    }
    if (class(mod)[1] == "phylolm") {
      if (!is.object(mod.r)) {
        Y <- mod$y
        mod.r <- lm(Y ~ 1)
      }
      if (!is.element(class(mod.r)[1], c("phylolm", "lm"))) {
        stop("mod.r must be class phylolm or lm.")
      }
      return(R2.lr.phylolm.original(mod, mod.r))
    }
    if (class(mod)[1] == "phyloglm") {
      if (!is.object(mod.r)) {
        Y <- mod$y
        mod.r <- glm(Y ~ 1, family = "binomial")
      }
      if (!is.element(class(mod.r)[1], c("phyloglm", "glm"))) {
        stop("mod.r must be class phyloglm or glm.")
      }
      return(R2.lr.phyloglm.original(mod, mod.r))
    }
  }
  
  R2.lr.lmerMod.original <- function(mod = NA, mod.r = NA) {
    
    X <- model.matrix(mod)
    n <- dim(X)[1]
    
    R2.lr <- 1 - exp(-2/n * (logLik(mod) - logLik(mod.r)))
    return(R2.lr)
  }
  
  R2.lr.glmerMod.original <- function(mod = NA, mod.r = NA) {
    
    X <- model.matrix(mod)
    n <- dim(X)[1]
    
    R2.lr <- (1 - exp(-2/n * (logLik(mod) - logLik(mod.r))))/(1 - exp(2/n * logLik(mod.r)))
    return(R2.lr)
  }
  
  R2.lr.phylolm.original <- function(mod = NA, mod.r = NA) {
    
    X <- mod$X
    n <- dim(X)[1]
    
    R2.lr <- 1 - exp(-2/n * (logLik(mod)[[1]] - logLik(mod.r)[[1]]))
    return(R2.lr)
  }
  
  R2.lr.phyloglm.original <- function(mod = NA, mod.r = NA) {
    
    Y <- mod$y
    n <- dim(mod$X)[1]
    X <- mod$X[, -1]
    
    alpha.cutoff <- 40
    if (mod$alpha < alpha.cutoff) {
      LL <- mod$logLik
    } else {
      LL <- logLik(glm(Y ~ X, family = "binomial"))
    }
    if (class(mod.r)[1] == "phyloglm") {
      if (mod.r$alpha < alpha.cutoff) {
        LL.r <- mod.r$logLik
      } else {
        X.r <- mod.r$X
        LL.r <- logLik(glm(Y ~ 0 + X.r, family = "binomial"))
      }
    } else {
      LL.r <- logLik(mod.r)
    }
    
    R2.lr <- (1 - exp(-2/n * (LL - LL.r)))/(1 - exp(2/n * LL.r))
    return(R2.lr[1])
  }
  
  #### end of original code ---------
  
  ## data ---
  # data 
  set.seed(123)
  p1 <- 10; nsample <- 10; n <- p1 * nsample
  d <- data.frame(x1 = rnorm(n = n), 
                  x2 = rnorm(n = n), 
                  u1 = rep(1:p1, each = nsample), 
                  u2 = rep(1:p1, times = nsample))
  d$u1 <- as.factor(d$u1); d$u2 <- as.factor(d$u2)
  
  # LMM: y with random intercept
  b1 <- 1; b2 <- -1; sd1 <- 1.5
  d$y_re_intercept <- b1 * d$x1 + b2 * d$x2 + 
    rep(rnorm(n = p1, sd = sd1), each = nsample) +  # random intercept u1
    rep(rnorm(n = p1, sd = sd1), times = nsample) + # random intercept u2
    rnorm(n = n)
  
  # LMM: y with random slope
  b1 <- 0; sd1 <- 1; sd.x1 <- 2
  d$y_re_slope <- b1 * d$x1 + 
    rep(rnorm(n = p1, sd = sd1), each = nsample) + # random intercept u1
    d$x1 * rep(rnorm(n = p1, sd = sd.x1), times = nsample) + # random slope u1
    rnorm(n = n)
  
  # GLMM
  b1 <- 1; sd1 <- 1.5
  prob <- rr2::inv.logit(b1 * d$x1 + rep(rnorm(n = p1, sd = sd1), each = nsample)) # random intercept u1
  d$y_binary <- rbinom(n = n, size = 1, prob = prob)
  
  # y for PGLS
  b1 <- 1.5; signal <- 0.7
  phy <- ape::compute.brlen(ape::rtree(n = n), method = "Grafen", power = 1)
  phy.x <- ape::compute.brlen(phy, method = "Grafen", power = .0001)
  x_trait <- ape::rTraitCont(phy.x, model = "BM", sigma = 1)
  e <- signal^0.5 * ape::rTraitCont(phy, model = "BM", sigma = 1) + (1-signal)^0.5 * rnorm(n=n)
  d$x_trait <- x_trait[match(names(e), names(x_trait))]
  d$y_pgls <- b1 * x_trait + e
  rownames(d) <- phy$tip.label	
  
  # y for Phylogenetic logistic regression
  b1 <- 1.5; signal <- 2
  e <- signal * ape::rTraitCont(phy, model = "BM", sigma = 1)
  e <- e[match(phy$tip.label, names(e))]
  d$y_phy_binary <- rbinom(n = n, size = 1, prob = rr2::inv.logit(b1 * d$x1 + e))
  
  # models to test
  
  # LMM ===
  z.f.lmm <- lme4::lmer(y_re_intercept ~ x1 + x2 + (1 | u1) + (1 | u2), data = d, REML = F)
  z.x.lmm <- lme4::lmer(y_re_intercept ~ x1 + (1 | u1) + (1 | u2), data = d, REML = F)
  z.v.lmm <- lme4::lmer(y_re_intercept ~ 1 + (1 | u2), data = d, REML = F)
  z.0.lmm <- lm(y_re_intercept ~ 1, data = d)
  
  expect_equal(R2.ls(z.f.lmm, z.x.lmm), R2.ls.original(z.f.lmm, z.x.lmm))
  expect_equal(rr2::R2.lr(z.f.lmm, z.x.lmm), R2.lr.original(z.f.lmm, z.x.lmm))
  expect_equal(R2.ce(z.f.lmm, z.x.lmm), R2.ce.original(z.f.lmm, z.x.lmm))

  expect_equal(R2.ls(z.f.lmm, z.v.lmm), R2.ls.original(z.f.lmm, z.v.lmm))
  expect_equal(R2.lr(z.f.lmm, z.v.lmm), R2.lr.original(z.f.lmm, z.v.lmm))
  expect_equal(R2.ce(z.f.lmm, z.v.lmm), R2.ce.original(z.f.lmm, z.v.lmm))

  expect_equal(R2.ls(z.f.lmm, z.0.lmm), R2.ls.original(z.f.lmm, z.0.lmm))
  expect_equal(R2.lr(z.f.lmm, z.0.lmm), R2.lr.original(z.f.lmm, z.0.lmm))
  expect_equal(R2.ce(z.f.lmm, z.0.lmm), R2.ce.original(z.f.lmm, z.0.lmm))

  expect_equal(R2.ls(z.f.lmm, z.0.lmm), R2.ls.original(z.f.lmm))  # by default
  expect_equal(R2.lr(z.f.lmm, z.0.lmm), R2.lr.original(z.f.lmm))
  expect_equal(R2.ce(z.f.lmm, z.0.lmm), R2.ce.original(z.f.lmm))

  # REML is updated to ML for R2.lr, but not R2.ls or R2.ce
  z.f.lmm.reml <- lme4::lmer(y_re_intercept ~ x1 + x2 + (1 | u1) + (1 | u2), data = d, REML = T)
  expect_warning(R2.lr(z.f.lmm.reml), "mod updated with REML=F")

  # Here is an example with a random slope
  z.f.lmm2 <- lme4::lmer(y_re_slope ~ 1 + (1 | u1) + (0 + x1 | u1), data = d, REML = F)
  z.s.lmm2 <- lme4::lmer(y_re_slope ~ 1 + (1 | u1), data = d, REML = F)

  expect_equal(R2.ls(z.f.lmm2, z.s.lmm2), R2.ls.original(z.f.lmm2, z.s.lmm2))
  expect_equal(R2.lr(z.f.lmm2, z.s.lmm2), R2.lr.original(z.f.lmm2, z.s.lmm2))
  expect_equal(R2.ce(z.f.lmm2, z.s.lmm2), R2.ce.original(z.f.lmm2, z.s.lmm2))

  # GLMM ===
  z.f.glmm <- lme4::glmer(y_binary ~ x1 + (1 | u1), data = d, family = "binomial")
  z.x.glmm <- lme4::glmer(y_binary ~ 1 + (1 | u1), data = d, family = "binomial")
  z.v.glmm <- glm(y_binary ~ x1, data = d, family = "binomial")

  expect_equal(R2.ls(z.f.glmm, z.x.glmm), R2.ls.original(z.f.glmm, z.x.glmm))
  expect_equal(R2.lr(z.f.glmm, z.x.glmm), R2.lr.original(z.f.glmm, z.x.glmm))
  expect_equal(R2.ce(z.f.glmm, z.x.glmm), R2.ce.original(z.f.glmm, z.x.glmm))

  expect_equal(R2.ls(z.f.glmm, z.v.glmm), R2.ls.original(z.f.glmm, z.v.glmm))
  expect_equal(R2.lr(z.f.glmm, z.v.glmm), R2.lr.original(z.f.glmm, z.v.glmm))
  expect_equal(R2.ce(z.f.glmm, z.v.glmm), R2.ce.original(z.f.glmm, z.v.glmm))

  expect_equal(R2.ls(z.f.glmm), R2.ls.original(z.f.glmm))  # by default
  expect_equal(R2.lr(z.f.glmm), R2.lr.original(z.f.glmm))
  expect_equal(R2.ce(z.f.glmm), R2.ce.original(z.f.glmm))

  # PGLS ===
  z.x.pgls <- phylolm::phylolm(y_pgls ~ 1, phy = phy, data = d, model = "lambda")
  lam.x.pgls <- round(z.x.pgls$optpar, digits = 4)
  z.f.pgls <- phylolm::phylolm(y_pgls ~ x_trait, phy = phy, data = d, model = "lambda", starting.value = 0.98 * lam.x.pgls + 0.01)
  z.v.pgls <- lm(y_pgls ~ x_trait, data = d)

  expect_equal(R2.ls(z.f.pgls, z.x.pgls, phy = phy), R2.ls.original(z.f.pgls, z.x.pgls, phy = phy))
  expect_equal(R2.lr(z.f.pgls, z.x.pgls), R2.lr.original(z.f.pgls, z.x.pgls))
  expect_equal(R2.ce(z.f.pgls, z.x.pgls, phy = phy), R2.ce.original(z.f.pgls, z.x.pgls, phy = phy))

  expect_equal(R2.ls(z.f.pgls, z.v.pgls, phy = phy), R2.ls.original(z.f.pgls, z.v.pgls, phy = phy))
  expect_equal(R2.lr(z.f.pgls, z.v.pgls), R2.lr.original(z.f.pgls, z.v.pgls))
  expect_equal(R2.ce(z.f.pgls, z.v.pgls, phy = phy), R2.ce.original(z.f.pgls, z.v.pgls, phy = phy))

  expect_equal(R2.ls(z.f.pgls, phy = phy), R2.ls.original(z.f.pgls, phy = phy))  # by default
  expect_equal(R2.lr(z.f.pgls), R2.lr.original(z.f.pgls))
  expect_equal(R2.ce(z.f.pgls, phy = phy), R2.ce.original(z.f.pgls, phy = phy))

  # P Log ===
  z.f.plog <- rr2::binaryPGLMM(y_phy_binary ~ x1, data = d, phy = phy)
  z.x.plog <- rr2::binaryPGLMM(y_phy_binary ~ 1, data = d, phy = phy)
  z.v.plog <- glm(y_phy_binary ~ x1, data = d, family = "binomial")

  expect_equal(R2.ls(z.f.plog, z.x.plog), R2.ls.original(z.f.plog, z.x.plog))
  expect_equal(R2.ce(z.f.plog, z.x.plog), R2.ce.original(z.f.plog, z.x.plog))

  expect_equal(R2.ls(z.f.plog, z.v.plog), R2.ls.original(z.f.plog, z.v.plog))
  expect_equal(R2.ce(z.f.plog, z.v.plog), R2.ce.original(z.f.plog, z.v.plog))

  expect_equal(R2.ls(z.f.plog), R2.ls.original(z.f.plog))  # by default
  expect_equal(R2.ce(z.f.plog), R2.ce.original(z.f.plog))

  # R.lr
  z.f.plog2 <- rr2::phyloglm(y_phy_binary ~ x1, data = d, start.alpha = 1, phy = phy, opt.method = "Nelder-Mead")
  z.x.plog2 <- rr2::phyloglm(y_phy_binary ~ 1, data = d, phy = phy, start.alpha = min(20, z.f.plog2$alpha), opt.method = "Nelder-Mead")
  z.v.plog2 <- glm(y_phy_binary ~ x1, data = d, family = "binomial")

  expect_equal(R2.lr(z.f.plog2, z.x.plog2), R2.lr.original(z.f.plog2, z.x.plog2))
  expect_equal(R2.lr(z.f.plog2, z.v.plog2), R2.lr.original(z.f.plog2, z.v.plog2))
  expect_equal(R2.lr(z.f.plog2), R2.lr.original(z.f.plog2))
})
