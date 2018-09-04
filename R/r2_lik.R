#' Calculate R2.lik
#'
#' Calculate partial and total R2s for LMM, GLMM, PGLS, and PGLMMs using R2.lik, an R2 based on the likelihood of observing the data.
#' 
#' @param mod A model with the following class: 'lmerMod', 'glmerMod', 'phylolm', 'phyloglm', or 'communityPGLMM'.
#' @param mod.r A reduced model. If a reduced model is not provided, the total R2 will be computed using the corresponding reduced model with the intercept as the only predictor.
#' @return R2.lik
#' @export
#'
#' @details  R2.lik works with classes lmerMod (LMM), glmerMod (GLMM), phylolm (PGLS), phyloglm (PGLMM), and communityPGLMM (Gaussian PLMM). 
#' 
#' \deqn{partial R2 = 1 - exp(-2/n*(logLik(mod.f) - logLik(mod.r)))}
#' 
#' where mod.f and mod.r are the full and reduced models, respectively. The total R2 is given when mod.r is the model corresponding to mod.f that contains only the intercept.
#' 
#' R2.lik is also computed for LMMs and GLMMs in the {MuMIn} package.
#' 
#' @seealso MuMIn
#' @author Anthony R. Ives
#' @references Ives A. in press. R2s for Correlated Data: Phylogenetic Models, LMMs, and GLMMs. Systematic Biology.
#' 
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
#' d$y <- b1 * d$x1 + b2 * d$x2 + rep(rnorm(n=p1, sd=sd1), each=nsample) + rep(rnorm(n=p1, sd=sd1), times=nsample) + rnorm(n=n)
#' 
#' z.f <- lmer(y ~ x1 + x2 + (1 | u1) + (1 | u2), data=d, REML = F)
#' z.x <- lmer(y ~ x1 + (1 | u1) + (1 | u2), data=d, REML = F)
#' z.v <- lmer(y ~ 1 + (1 | u2), data=d, REML = F)
#' z.0 <- lm(y ~ 1, data=d)
#' 
#' R2.lik(z.f, z.x)
#' R2.lik(z.f, z.v)
#' R2.lik(z.f)
#' 
#' # These give the same results
#' R2.lik(z.f, z.0)
#' R2.lik(z.f)
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
#' R2.lik(z.f, z.x)
#' R2.lik(z.f, z.v)
#' R2.lik(z.f)
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
#' R2.lik(z.f, z.x)
#' R2.lik(z.f, z.v)
#' R2.lik(z.f)
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
#' z.f <- phyloglm(y ~ x, data=d, start.alpha = 1, phy=phy)
#' z.x <- phyloglm(y ~ 1, data=d, phy=phy, start.alpha=min(20,z.f$alpha))
#' z.v <- glm(y ~ x, data=d, family="binomial")
#' 
#' R2.lik(z.f, z.x)
#' R2.lik(z.f, z.v)
#' R2.lik(z.f)
#' 


R2.lik <- function(mod = NULL, mod.r = NULL) {
  if(class(mod)[1] == "merModLmerTest") 
    class(mod) = "lmerMod"
  
  if (!is.element(class(mod)[1], c("lm", "glm", "lmerMod", "glmerMod", "phylolm", "phyloglm", "communityPGLMM"))) {
    stop("mod must be class one of classes lm, glm, lmerMod, glmerMod, phylolm, phyloglm, communityPGLMM.")
  }
  
  if (class(mod)[1] == "lm") {
    if (!is.object(mod.r)) {
      Y <- model.frame(mod)[, 1]
      mod.r <- lm(Y ~ 1)
    }
    if (!is.element(class(mod.r)[1], c("lm"))) {
      stop("mod.r must be class lm.")
    }
    return(R2.lik.lm(mod, mod.r)[1])
  }
  
  if (class(mod)[1] == "glm") {
    if (!is.object(mod.r)) {
      Y <- model.frame(mod)[, 1]
      mod.r <- glm(Y ~ 1, family = family(mod)[[1]])
    }
    if (!is.element(class(mod.r)[1], c("glm"))) {
      stop("mod.r must be class glm.")
    }
    return(R2.lik.glm(mod, mod.r)[1])
  }
  
  if (class(mod)[1] == "lmerMod") {
    if (!is.object(mod.r)) {
      Y <- model.frame(mod)[, 1]
      mod.r <- lm(Y ~ 1)
    }
    if(class(mod.r)[1] == "merModLmerTest") class(mod.r) = "lmerMod"
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
    return(R2.lik.lmerMod(mod, mod.r)[1])
  }
  
  if (class(mod)[1] == "glmerMod") {
    if (!is.object(mod.r)) {
      Y <- model.frame(mod)[, 1]
      mod.r <- glm(Y ~ 1, family = family(mod)[[1]])
    }
    if (!is.element(class(mod.r)[1], c("glmerMod", "glm"))) {
      stop("mod.r must be class glmerMod or glm.")
    }
    return(R2.lik.glmerMod(mod, mod.r)[1])
  }
  
  if (class(mod)[1] == "phylolm") {
    if (!is.object(mod.r)) {
      Y <- mod$y
      mod.r <- lm(Y ~ 1)
    }
    if (!is.element(class(mod.r)[1], c("phylolm", "lm"))) {
      stop("mod.r must be class phylolm or lm.")
    }
    return(R2.lik.phylolm(mod, mod.r))
  }
  
  if (class(mod)[1] == "phyloglm") {
    if (!is.object(mod.r)) {
      Y <- mod$y
      mod.r <- glm(Y ~ 1, family = "binomial")
    }
    if (!is.element(class(mod.r)[1], c("phyloglm", "glm"))) {
      stop("mod.r must be class phyloglm or glm.")
    }
    return(R2.lik.phyloglm(mod, mod.r)[1])
  }
  
  if (class(mod)[1] == "communityPGLMM") {
    if(mod$family == "binomial")
      stop("binary communityPGLMMs do not have log likelihood, \
           If you are interested in LRT of random terms, use \
           phyr::communityPGLMM.binary.LRT()")
    if(mod$REML == TRUE)
      stop("mod was fitted with REML, please set it to FALSE and re-fit it")
    
    if (!is.object(mod.r)) {
      Y <- mod$Y
      mod.r <- lm(Y ~ 1)
    }
    if (!is.element(class(mod.r)[1], c("communityPGLMM", "lm"))) {
      stop("mod.r must be class communityPGLMM or lm.")
    }
    return(R2.like.communityPGLMM(mod, mod.r))
  }
}

R2.lik.lm <- function(mod = NULL, mod.r = NULL) {
  X <- model.matrix(mod)
  n <- dim(X)[1]
  
  R2.lik <- 1 - exp(-2/n * (logLik(mod)[[1]] - logLik(mod.r)[[1]]))
  return(R2.lik)
}

R2.lik.glm <- R2.lik.lm

R2.lik.lmerMod <- R2.lik.lm

R2.lik.glmerMod <- function(mod = NULL, mod.r = NULL) {
  X <- model.matrix(mod)
  n <- dim(X)[1]
  R2.lik <- (1 - exp(-2/n * (logLik(mod) - logLik(mod.r))))/(1 - exp(2/n * logLik(mod.r)))
  return(R2.lik)
}

R2.lik.phylolm <- function(mod = NULL, mod.r = NULL) {
  X <- mod$X
  n <- dim(X)[1]
  R2.lik <- 1 - exp(-2/n * (logLik(mod)[[1]] - logLik(mod.r)[[1]]))
  return(R2.lik)
}

R2.lik.phyloglm <- function(mod = NULL, mod.r = NULL) {
  
  Y <- mod$y
  X <- mod$X
  n <- dim(X)[1]
  
  alpha.cutoff <- 40
  if (mod$alpha < alpha.cutoff) {
    LL <- mod$logLik
  } else {
    LL <- logLik(glm(Y ~ 0 + X, family = "binomial"))
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
  
  R2.lik <- (1 - exp(-2/n * (LL - LL.r)))/(1 - exp(2/n * LL.r))
  
  return(R2.lik)
}

R2.like.communityPGLMM <- function(mod = NULL, mod.r = NULL){
  n <- nrow(mod$X)
  if(class(mod.r) == "lm"){
    ll.r <- logLik(mod.r)[[1]]
  } else {
    ll.r <- mod.r$logLik
  }
  R2.lik <- 1 - exp(-2/n * (mod$logLik - ll.r))
  return(R2.lik)
}
