#' Calculate R2.resid
#'
#' Calculate partial and total R2s for LMM, GLMM, PGLS, and PGLMM using R2.resid, an extension of ordinary least-squares (OLS) R2s. For LMMs and GLMMs, R2.resid is related to the method proposed by Nakagawa and Schielzeth (2013).
#' 

#' @param mod A regression model with one of the following classes: 'lm', 'glm', 'lmerMod', 'glmerMod', 'phylolm', 'gls', or 'binaryPGLMM'. For 'glmerMod', only family = c('binomial', 'poisson') are supported.
#' @param mod.r A reduced model; if not provided, the total R2 will be given by setting 'mod.r' to the model corresponding to 'mod' with the intercept as the only predictor.
#' @param phy The phylogeny for phylogenetic models (as a 'phylo' object), which must be specified for models of class `phylolm`.
#' @param sigma2_d Distribution-specific variance \eqn{\sigma^2_d}{sigma2d} (see Details). For binomial GLMs, GLMMs and PGLMMs with logit link functions, options are c('s2w', 'NS', 'rNS'). For binomial GLMs, GLMMs and PGLMMs with probit link functions, options are c('s2w', 'NS'). Other families use 's2w'.
#' 
#' @details  R2.resid works with classes 'lm', 'glm', 'lmerMod', 'glmerMod', 'phylolm', and 'binaryPGLMM'.
#' 
#' \strong{LMM (lmerMod):}
#' 
#' \deqn{partial R^2 = 1 - \sigma^2_{e.f}/\sigma^2_{e.r}}{partial R2 = 1 - σ2e.f/σ2e.r}
#'    
#' \deqn{total R^2 = 1 - \sigma^2_{e.f}/var(y)}{total R2 = 1 - σ2e.f/var(y)}
#'    
#' where \eqn{\sigma^2_{e.f}}{σ2e.f} and \eqn{\sigma^2_{e.r}}{σ2e.r} are the estimated residual variances from the full and reduced LMM, and var(y) is the total variance of the response (dependent) variable.
#' 
#' \strong{GLMM (glmerMod):} 
#' 
#' \deqn{total R^2 = 1 - \sigma^2_d/(\sigma^2_x + \sigma^2_b + \sigma^2_d)}{total R2 = 1 - σ2d/(σ2x + σ2b + σ2d)}
#'    
#' where \eqn{\sigma^2_x}{σ2x} and \eqn{\sigma^2_b}{σ2b} are the estimated variances associated with the fixed and random effects. \eqn{\sigma^2_d}{σ2d} is a term that scales the implied 'residual variance' of the GLMM (see Ives 2018, Appendix 1). The default used for \eqn{\sigma^2_d}{σ2d} is \eqn{\sigma^2_w}{σ2w} which is computed from the iterative weights of the GLMM. Specifically,
#' 
#' \deqn{\sigma_{w}^{2}=var(g'(\mu)*(y-\mu))}{σ2w = var(g'(μ) * (y – μ))}
#' 
#' where g'() is the derivative of the link function, and \eqn{(y-\mu)}{(y – μ)} is the difference between the data y and their predicted values \eqn{\mu}{μ}. This is the default option specified by sigma2_d = 's2w'. For binomial models with a logit link function, sigma2_d = 'NS' gives the scaling \eqn{\sigma^2_d =  \pi^2/3}{σ2d =  π^2/3}, and sigma2_d = 'rNS' gives \eqn{\sigma^2_d = 0.8768809 * \pi^2/3}{σ2d = 0.8768809 * π^2/3}. For binomial models with a probit link function, sigma2_d = 'NS' gives the scaling \eqn{\sigma^2_d = 1}{σ2d = 1}. In general option sigma2_d = 's2w' will give values lower than sigma2_d = 'NS' and 'rNS', but the values will be closer to \code{R2.lik()} and \code{R2.pred()}. For other forms of sigma2_d from Nakagawa and Schielzeth (2013) and Nakagawa et al. (2017), see the MuMIn package.
#' 
#' Partial R2s are given by the standard formula
#' 
#' \deqn{partial R^2 = 1 - (1 - R^2_{.f})/(1 - R^2_{.r})}{partial R2 = 1 - (1 - R2.f)/(1 - R2.r)}
#' 
#' where R2.f and R2.r are the total R2s for full and reduced models, respectively.
#' 
#' \strong{PGLS (phylolm):} 
#' 
#' \deqn{partial R^2 = 1 - c.f * \sigma^2_{.f}/(c.r * \sigma^2_{.r})}{partial R2 = 1 - c.f * σ2.f/(c.r * σ2.r)}
#' 
#' where \eqn{\sigma^2_{.f}}{σ2.f} and \eqn{\sigma^2_{.r}}{σ2.r} are the variances estimated for the PGLS full and reduced models, and c.f and c.r are the scaling values for full and reduce models that equal the total sum of phylogenetic branch length estimates. Note that the phylogeny needs to be specified in R2.resid.
#' 
#' Note that \code{phylolm()} can have difficulties in finding solutions when there is no phylogenetic signal;
#' when the estimate indicates no phylogenetic signal, you should refit the model with the corresponding LM.
#' 
#' \strong{PGLMM (binaryPGLMM):} 
#' 
#' The binary PGLMM is computed in the same way as the binomial GLMM, with options sigma_d = c('s2w', 'NS', 'rNS'). The estimated variance of the random effect associated with the phylogeny, \eqn{\sigma^2_b}{σ2b}, is multiplied by the diagonal elements of the phylogenetic covariance matrix. For binary models, this covariance matrix should be standardized so that all diagonal elements are the same (a contemporaneous or ultrametric phylogenetic tree) (Ives and Garland 2014). In case this is not done, however, the code takes the geometric average of the diagonal elements.
#' 
#' Note that the version of \code{binaryPGLMM()} in the package ape is replaced by a version contained within {rr2} that outputs all of the required information for the calculation of \code{R2.resid()}
#' 
#' \strong{LM (lm) and GLM (glm):} 
#' 
#' For compatibility and generating reduced models, rr2 will compute \code{R2.resid()} for LM and GLM that correspond to LMM/PGLS and GLMM/PGLMM.
#' 
#' @return R2.resid value.
#' @importFrom lme4 VarCorr
#' @seealso MuMIn, lme4, ape, phylolm, pez
#' @export
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
#' R2.resid(z.f, z.x)
#' R2.resid(z.f, z.v)
#' R2.resid(z.f)
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
#' R2.resid(z.f, z.x)
#' R2.resid(z.f, z.v)
#' R2.resid(z.f)
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
#' R2.resid(z.f, z.x, phy = phy)
#' R2.resid(z.f, z.v, phy = phy)
#' R2.resid(z.f, phy = phy)
#' 
#' # This also works for models fit with gls() in {nlme}
#' z.x <- gls(y ~ 1, data = d, correlation = corPagel(1, phy), method = "ML")
#' z.f <- gls(y ~ x, data = d, correlation = corPagel(1, phy), method = "ML")
#' z.v <- lm(y ~ x, data = d)
#' 
#' R2.resid(z.f, z.x)
#' R2.resid(z.f, z.v)
#' R2.resid(z.f)
#' 
#' # But note that you need to define weights for gls() with non-ultrametric trees; 
#' # if not, you will get a error "Matrix is not block-diagonal"
#' 
#' phy.nu <- rtree(n = n)
#' 
#' # Generate random data
#' e <- signal^0.5 * rTraitCont(phy.nu, model = 'BM', sigma = 1) + (1-signal)^0.5 * rnorm(n = n)
#' d$x <- x[match(names(e), names(x))]
#' d$y <- b1 * x + e
#' rownames(d) <- phy.nu$tip.label
#' 
#' weights <- diag(vcv.phylo(phy.nu))
#' z.x <- gls(y ~ 1,data = d, 
#'          correlation = corPagel(1, phy.nu),
#'          weights=varFixed(~weights), method = "ML")
#' z.f <- gls(y ~ x,data = d, 
#'          correlation = corPagel(1, phy.nu),
#'          weights=varFixed(~weights), method = "ML")
#' z.v <- lm(y ~ x, data = d)
#' 
#' R2.resid(z.f, z.x)
#' R2.resid(z.f, z.v)
#' R2.resid(z.f)
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
#' R2.resid(z.f, z.x, phy = phy)
#' R2.resid(z.f, z.v, phy = phy)
#' R2.resid(z.f, phy = phy)
#' 
#' @author Anthony R. Ives
#' @seealso MuMIn
#' @references Ives A.R. and Li D. 2018. rr2: An R package to calculate R2s for regression models. Journal of Open Source Software. DOI:10.21105/joss.01028
#' 
#' Ives A.R. 2018. R2s for Correlated Data: Phylogenetic Models, LMMs, and GLMMs. Systematic Biology. DOI:10.1093/sysbio/syy060
#' 
#' Ives A. R., Garland T., Jr. 2014. Phylogenetic regression for binary dependent variables. In: Garamszegi LZ editor. Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology. Berlin Heidelberg, Springer-Verlag, p. 231-261.
#' 
#' Nakagawa S., Schielzeth H. 2013. A general and simple method for obtaining R2 from generalized linear mixed-effects models. Methods in Ecology and Evolution, 4:133-142.
#' 
#' Nakagawa S., Johnson P. C. D., Schielzeth H. 2017. The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. Journal of the Royal Society Interface, 14.
#' 
R2.resid <- function(mod = NULL, mod.r = NULL, phy = NULL, 
                     sigma2_d = c("s2w", "NS", "rNS")) {
    if (class(mod)[1] == "merModLmerTest") 
        class(mod) <- "lmerMod"
    
    if (!is.element(class(mod)[1], c("lm", "glm", "lmerMod", "glmerMod", "phylolm", "gls", "binaryPGLMM"))) {
        stop("mod must be class one of classes lm, glm, lmerMod, glmerMod, phylolm, gls, binaryPGLMM.")
    }
    
    sigma2_d <- match.arg(sigma2_d)
    if (!is.null(sigma2_d) && !is.element(sigma2_d, c("s2w", "NS", "rNS"))) 
        stop("Please specify residual variance c('s2w', 'NS', 'rNS').")
    if (is.null(sigma2_d)) 
        sigma2_d <- "s2w"
    
    if (class(mod)[1] == "lm") {
        if (!is.object(mod.r)) {
            y <- model.frame(mod)[, 1]
            mod.r <- lm(y ~ 1)
        }
        if (!is.element(class(mod.r)[1], c("lm"))) {
            stop("mod.r must be class lm.")
        }
        return(R2.resid.lm(mod, mod.r))
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
        if (!is.element(family(mod)[[1]], c("binomial", "poisson"))) {
            stop("Sorry, but R2.resid only works for family = binomial or poisson.")
        }
        
        return(R2.resid.glm(mod, mod.r, sigma2_d = sigma2_d))
    }
    
    if (class(mod)[1] == "lmerMod") {
        if (!is.object(mod.r)) {
            # exists()?
            y <- model.frame(mod)[, 1]
            mod.r <- lm(y ~ 1)
        }
        if (class(mod.r)[1] == "merModLmerTest") 
            class(mod.r) <- "lmerMod"
        if (!is.element(class(mod)[1], c("lmerMod", "lm"))) {
            stop("mod.r must be class lmerMod or lm.")
        }
        return(R2.resid.lmerMod(mod, mod.r))
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
        if (!is.element(family(mod)[[1]], c("binomial", "poisson"))) {
            stop("Sorry, but R2.resid only works for family = binomial or poisson.")
        }
        
        return(R2.resid.glmerMod(mod, mod.r, sigma2_d = sigma2_d))
    }
    
    if (class(mod)[1] == "phylolm") {
      if (!is.object(phy)) {
        stop("For phylolm you must provide the phylo object")
      }
      if (!is.object(mod.r)) {
        y <- mod$y
        mod.r <- lm(y ~ 1)
      }
      if (!is.element(class(mod.r)[1], c("phylolm", "lm"))) {
        stop("mod.r must be class phylolm or lm.")
      }
      return(R2.resid.phylolm(mod, mod.r, phy))
    }
    
    if (class(mod)[1] == "gls") {
      if (!is.object(mod.r)) {
        y <- as.numeric(fitted(mod)+resid(mod))
        mod.r <- lm(y ~ 1)
      }
      if (!is.element(class(mod.r)[1], c("gls", "lm"))) {
        stop("mod.r must be class gls or lm.")
      }
      return(R2.resid.gls(mod, mod.r))
    }
    
    if (class(mod)[1] == "binaryPGLMM") {
        if (!is.object(mod.r)) {
            y <- mod$y
            mod.r <- glm(y ~ 1, family = "binomial")
        }
        if (!is.element(class(mod.r)[1], c("binaryPGLMM", "glm"))) {
            stop("mod.r must be class binaryPGLMM or glm.")
        }
        
        return(R2.resid.binaryPGLMM(mod, mod.r, sigma2_d = sigma2_d))
    }
}

R2.resid.lm <- function(mod = NULL, mod.r = NULL) {
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

R2.resid.glm <- function(mod = NULL, mod.r = NULL, sigma2_d = sigma2_d) {
    mu <- mod$fitted.values
    Yhat <- family(mod)$linkfun(mu)
    y <- model.frame(mod)[,1]
    if (family(mod)[1] == "binomial") {
        if (is.matrix(model.frame(mod)[, 1])) {
            size <- rowSums(model.frame(mod)[, 1])
            y <- y[,1]/size
        } else {
            size <- 1
        }
        if (family(mod)[2] == "logit") {
            if (sigma2_d == "s2w") 
              sig2e <- var((y - mu)/(mu * (1 - mu)))
            if (sigma2_d == "NS") 
                sig2e <- pi^2/3
            if (sigma2_d == "rNS") 
                sig2e <- 0.8768809 * pi^2/3
        } else {
            if (sigma2_d == "s2w") 
                sig2e <- var((y - mu)/dnorm(qnorm(mu)))
            if (sigma2_d == "NS") 
                sig2e <- 1
            if (sigma2_d == "rNS") 
                sig2e <- 1
        }
    }
    if (family(mod)[1] == "poisson") 
      sig2e <- var((y - mu)/mu)
    
    SSE.resid <- sig2e/(var(Yhat) + sig2e)
    
    mu.r <- mod.r$fitted.values
    Yhat.r <- family(mod.r)$linkfun(mu.r)
    if (family(mod.r)[1] == "binomial") 
        if (family(mod.r)[2] == "logit") {
            if (sigma2_d == "s2w") 
                sig2e.r <- var(size * (y - mu.r)/(mu.r * (1 - mu.r)))
            if (sigma2_d == "NS") 
                sig2e.r <- pi^2/3
            if (sigma2_d == "rNS") 
                sig2e.r <- 0.8768809 * pi^2/3
        } else {
            if (sigma2_d == "s2w") 
                sig2e.r <- var(size * (y - mu.r)/dnorm(qnorm(mu.r)))
            if (sigma2_d == "NS") 
                sig2e.r <- 1
            if (sigma2_d == "rNS") 
                sig2e.r <- 1
        }
    if (family(mod.r)[1] == "poisson") 
      sig2e.r <- var((y - mu.r)/mu.r)
    
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

R2.resid.glmerMod <- function(mod = NULL, mod.r = NULL, sigma2_d = sigma2_d) {
    
    X <- model.matrix(mod)
    y <- model.frame(mod)[, 1]
    n <- dim(X)[1]
    X.r <- model.matrix(mod.r)

    # full model
    mu <- fitted(mod)
    Yhat <- X %*% lme4::fixef(mod)
    if (family(mod)[1] == "binomial") {
        if (is.matrix(model.frame(mod)[, 1])) {
            size <- rowSums(model.frame(mod)[, 1])
            y <- y[,1]/size
        } else {
            size <- 1
        }
        if (family(mod)[2] == "logit") {
            if (sigma2_d == "s2w") 
              sig2e <- var((y - mu)/(mu * (1 - mu)))
            if (sigma2_d == "NS") 
                sig2e <- pi^2/3
            if (sigma2_d == "rNS") 
                sig2e <- 0.8768809 * pi^2/3
        } else {
            if (sigma2_d == "s2w") 
                sig2e <- var((y - mu)/dnorm(qnorm(mu)))
            if (sigma2_d == "NS") 
                sig2e <- 1
            if (sigma2_d == "rNS") 
                sig2e <- 1
        }
    }
    if (family(mod)[1] == "poisson") 
      sig2e <- var((y - mu)/mu)
    
    sig2a <- VarCorr(mod)[[1]][1]
    nranef <- length(VarCorr(mod))
    if (nranef > 1) 
        for (i in 2:nranef) sig2a <- sig2a + VarCorr(mod)[[i]][1]
    
    SSE.resid <- sig2e/(var(Yhat) + sig2a + sig2e)
    
    # reduced model
    if (class(mod.r)[1] == "glmerMod") {
        X.r <- model.matrix(mod.r)
        mu.r <- fitted(mod.r)
        Yhat.r <- X.r %*% lme4::fixef(mod.r)
        if (family(mod.r)[1] == "binomial") 
            if (family(mod.r)[2] == "logit") {
                if (sigma2_d == "s2w") 
                   sig2e.r <- var((y - mu.r)/(mu.r * (1 - mu.r)))
                if (sigma2_d == "NS") 
                  sig2e.r <- pi^2/3
                if (sigma2_d == "rNS") 
                  sig2e.r <- 0.8768809 * pi^2/3
            } else {
                if (sigma2_d == "s2w") 
                  sig2e.r <- var((y - mu.r)/dnorm(qnorm(mu.r)))
                if (sigma2_d == "NS") 
                  sig2e.r <- 1
                if (sigma2_d == "rNS") 
                  sig2e.r <- 1
            }
        if (family(mod.r)[1] == "poisson") 
           sig2e.r <- var((y - mu.r)/mu.r)
        
        sig2a.r <- VarCorr(mod.r)[[1]][1]
        nranef.r <- length(VarCorr(mod.r))
        if (nranef.r > 1) 
            for (i in 2:nranef.r) sig2a.r <- sig2a.r + VarCorr(mod.r)[[i]][1]
        
        SSE.resid.r <- sig2e.r/(var(Yhat.r) + sig2a.r + sig2e.r)
    }
    
    if (class(mod.r)[1] == "glm") {
      mu.r <- mod.r$fitted.values
        Yhat.r <- family(mod.r)$linkfun(mu.r)
        if (family(mod.r)[1] == "binomial") 
            if (family(mod.r)[2] == "logit") {
                if (sigma2_d == "s2w") 
                  sig2e.r <- var((y - mu.r)/(mu.r * (1 - mu.r)))
                if (sigma2_d == "NS") 
                  sig2e.r <- pi^2/3
                if (sigma2_d == "rNS") 
                  sig2e.r <- 0.8768809 * pi^2/3
            } else {
                if (sigma2_d == "s2w") 
                  sig2e.r <- var((y - mu.r)/dnorm(qnorm(mu.r)))
                if (sigma2_d == "NS") 
                  sig2e.r <- 1
                if (sigma2_d == "rNS") 
                  sig2e.r <- 1
            }
        if (family(mod.r)[1] == "poisson") 
          sig2e.r <- var((y - mu.r)/mu.r)
        
        SSE.resid.r <- sig2e.r/(var(Yhat.r) + sig2e.r)
    }
    
    R2.resid <- 1 - SSE.resid/SSE.resid.r
    return(R2.resid[1])
}

R2.resid.phylolm <- function(mod = NULL, mod.r = NULL, phy = NULL) {
    
    X <- mod$X
    n <- dim(X)[1]
    
    if (!mod$model %in% c("lambda", "OUrandomRoot", "OUfixedRoot", "BM", "kappa", 
        "delta", "EB", "trend")) {
        stop("Evolution model not supported yet.")
    }
    
    phy.f <- transf_phy(mod, phy)
    
    scal <- sum(phy.f$edge.length)/n
    sigma2 <- mod$sigma2
    
    if (class(mod.r) == "phylolm") {
        if (!mod.r$model %in% c("lambda", "OUrandomRoot", "OUfixedRoot", "BM", "kappa", 
            "delta", "EB", "trend")) {
            stop("Evolution model not supported yet.")
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

R2.resid.gls <- function(mod = NULL, mod.r = NULL) {
  
  n <- mod$dims$N
  
  cormatrix <- nlme::corMatrix(mod$modelStruct$corStruct)
  if(!is.null(attr(mod$modelStruct$varStruct, 'weights'))){
    VCVdiag <- 1/attr(mod$modelStruct$varStruct, 'weights')
    VCV.f <- diag(VCVdiag) %*% cormatrix %*% diag(VCVdiag)
    phy.f <- ape::vcv2phylo(VCV.f)
  }else{
    phy.f <- ape::vcv2phylo(cormatrix)
  }
  
  scal <- sum(phy.f$edge.length)/n
  sigma2 <- mod$sigma^2
  
  if (class(mod.r) == "gls") {

    cormatrix.r <- nlme::corMatrix(mod.r$modelStruct$corStruct)
    if(!is.null(attr(mod.r$modelStruct$varStruct, 'weights'))){
      VCVdiag <- 1/attr(mod.r$modelStruct$varStruct, 'weights')
      VCV.r <- diag(VCVdiag) %*% cormatrix.r %*% diag(VCVdiag)
      phy.r <- ape::vcv2phylo(VCV.r)
    }else{
      phy.r <- ape::vcv2phylo(cormatrix.r)
    }
    
    scal.r <- sum(phy.r$edge.length)/n
    sigma2.r <- mod.r$sigma^2
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

R2.resid.binaryPGLMM <- function(mod = NULL, mod.r = NULL, sigma2_d = sigma2_d) {
    
    y <- mod$y
    n <- length(y)
    Yhat <- mod$X %*% mod$B
    phyV <- mod$VCV
    s2 <- mod$s2
    scal <- prod(diag(s2 * phyV))^(1/n)
    mu <- mod$mu
    Yhat <- log(mu/(1 - mu))
    if (sigma2_d == "s2w") 
        sig2e <- var((y - mu)/(mu * (1 - mu)))
    if (sigma2_d == "NS") 
        sig2e <- pi^2/3
    if (sigma2_d == "rNS") 
        sig2e <- 0.8768809 * pi^2/3
    
    SSE.resid <- sig2e/(var(Yhat) + scal + sig2e)
    
    # reduced model
    if (class(mod.r)[1] == "binaryPGLMM") {
        Yhat.r <- mod.r$X %*% mod.r$B
        phyV.r <- mod.r$VCV
        s2.r <- mod.r$s2
        scal.r <- prod(diag(s2.r * phyV.r))^(1/n)
        mu.r <- mod.r$mu
        Yhat.r <- log(mu.r/(1 - mu.r))
        if (sigma2_d == "s2w") 
          sig2e.r <- var((y - mu.r)/(mu.r * (1 - mu.r)))
        if (sigma2_d == "NS") 
            sig2e.r <- pi^2/3
        if (sigma2_d == "rNS") 
            sig2e.r <- 0.8768809 * pi^2/3
        
        SSE.resid.r <- sig2e.r/(var(Yhat.r) + scal.r + sig2e.r)
    }
    
    if (class(mod.r)[1] == "glm") {
        mu.r <- mod.r$fitted.values
        Yhat.r <- log(mu.r/(1 - mu.r))
        if (sigma2_d == "s2w") 
            sig2e.r <- var((y - mu.r)/(mu.r * (1 - mu.r)))
        if (sigma2_d == "NS") 
            sig2e.r <- pi^2/3
        if (sigma2_d == "rNS") 
            sig2e.r <- 0.8768809 * pi^2/3
        
        SSE.resid.r <- sig2e.r/(var(Yhat.r) + sig2e.r)
    }
    
    R2.resid <- 1 - SSE.resid/SSE.resid.r
    
    return(R2.resid[1])
}
