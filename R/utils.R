#' @importFrom utils tail
#' @importFrom stats anova binomial family glm lm logLik model.frame model.matrix model.response optim pchisq pnorm poisson quantile reorder rnorm sd sigma update var na.omit fitted predict
#' @importMethodsFrom Matrix t %*% crossprod diag tcrossprod solve determinant update
NULL

#' Invert logit function
#'
#' Convert numeric values between 0 and 1.
#'
#' @param x a vector of numeric values
#' @export
inv.logit <- function(x) {
  1/(1 + exp(-x))
}

#' Partial R2
#'
#' Get partial R2 by comparing a model and its reduced model.
#'
#' @param mod a linear regression model
#' @param mod.r a reduced model based on \code{mod}
#' @return R2 value between 0 and 1
#' @export
partialR2 <- function(mod, mod.r) {
  anova.full <- anova(mod)
  anova.reduced <- anova(mod.r)
  
  sse.full <- tail(anova.full$`Sum Sq`, 1)
  sse.reduced <- tail(anova.reduced$`Sum Sq`, 1)
  
  return((sse.reduced - sse.full)/sse.reduced)
}

#' Adjusted partial R2
#'
#' Get adjusted partial R2 by comparing a model and its reduced model.
#'
#' @param mod a linear regression model
#' @param df.f degree of freedom of the \code{mod}
#' @param mod.r a reduced model based on \code{mod}
#' @param df.r degree of freedom of the reduced \code{mod.r}
#' @return A list of both R2 and adjusted R2 , the latter is not necessary to be between 0 and 1.
#' @export
partialR2adj <- function(mod, df.f, mod.r, df.r) {
  anova.full <- anova(mod)
  anova.reduced <- anova(mod.r)
  
  sse.full <- tail(anova.full$`Sum Sq`, 1)
  sse.reduced <- tail(anova.reduced$`Sum Sq`, 1)
  
  R2 <- 1 - sse.full/sse.reduced
  
  R2.adj <- 1 - (sse.full/df.f)/(sse.reduced/df.r)
  
  return(list(R2 = R2, R2.adj = R2.adj))
}

#' Transform a phylogeny based on a phylolm model
#' 
#' Using a fitted phylolm model to transform branch lengths of a phylogeny
#' 
#' @param phylolmMod a fitted phylolm model
#' @param phy a phylogeny with class "phylo"
#' @return a transformed phylogeny
#' @export
#' 
transf_phy <- function(phylolmMod, phy){
  if (!phylolmMod$model %in% c("BM", "trend")) {
    # optpar for BM models is NULL
    optpar <- round(phylolmMod$optpar, digits = 4)
    m.list <- list(x = optpar)
    
    if (phylolmMod$model %in% c("OUrandomRoot", "OUfixedRoot")) {
      names(m.list) <- "alpha"
    } else {
      names(m.list) <- phylolmMod$model
    }
    
    phy.f <- phylolm::transf.branch.lengths(phy, parameters = m.list, model = phylolmMod$model)$tree
  } else {
    # If model='BM' or model='trend', the output tree is the same as the input 
    # tree except that the output tree is in pruningwise order.
    phy.f <- phylolm::transf.branch.lengths(phy, parameters = NULL, model = phylolmMod$model)$tree
  }
  
  phy.f
}

#' Phylogenetic GLM for binary data
#'
#' Fitting phylogenetic generalized linear models for binary data (0 and 1).
#'
#' @param formula regression formula
#' @param data data frame to fit the model with
#' @param phy phylogenetic tree of type phylo with branch lengths
#' @param s2.init initial variance values for random terms, default is 0.1
#' @param B.init initial coefficient values for fixed terms, if not provided, will use those from \code{lm}
#' @param tol.pql tolerance value, default is 10^-6
#' @param maxit.pql the number of iterations, default is 200
#' @param maxit.reml the number of iterations for optim, default is 100
#' @return a large list with class as \code{binaryPGLMM}
#' @export
#'
binaryPGLMM <- function(formula, data = list(), phy, s2.init = 0.1, B.init = NULL, tol.pql = 10^-6, maxit.pql = 200, maxit.reml = 100) {
  
  # Helper function for \code{binaryPGLMM}
  
  # par = s2, tinvW = invW, tH = H, tVphy = Vphy, tX = X save(s2, invW, H, Vphy, X, file = 'pglmm.reml.RData')
  
  pglmm.reml <- function(par, tinvW, tH, tVphy, tX) {
    n <- dim(tX)[1]
    p <- dim(tX)[2]
    ss2 <- abs(Re(par))
    Cd <- ss2 * tVphy
    V <- tinvW + Cd
    LL <- 10^10
    if (sum(is.infinite(V)) == 0) {
      if (all(eigen(V)$values > 0)) {
        invV <- solve(V)
        logdetV <- determinant(V)$modulus[1]
        if (is.infinite(logdetV)) {
          cholV <- chol(V)
          logdetV <- 2 * sum(log(diag(chol(V))))
        }
        LL <- logdetV + t(tH) %*% invV %*% tH + determinant(t(tX) %*% invV %*% tX)$modulus[1]
      }
    }
    return(LL)
  }
  
  
  if (!inherits(phy, "phylo")) 
    stop("Object \"phy\" is not of class \"phylo\".")
  if (is.null(phy$edge.length)) 
    stop("The tree has no branch lengths.")
  if (is.null(phy$tip.label)) 
    stop("The tree has no tip labels.")
  phy <- reorder(phy, "postorder")
  n <- length(phy$tip.label)
  mf <- model.frame(formula = formula, data = data)
  if (nrow(mf) != length(phy$tip.label)) 
    stop("Number of rows of the design matrix does not match with length of the tree.")
  if (is.null(rownames(mf))) {
    warning("No tip labels, order assumed to be the same as in the tree.\n")
    data.names <- phy$tip.label
  } else {
    data.names <- rownames(mf)
  }
  .order <- match(data.names, phy$tip.label)  # do not name an object as a base function
  if (sum(is.na(.order)) > 0) {
    warning("Data names do not match with the tip labels.\n")
    rownames(mf) <- data.names
  } else {
    tmp <- mf
    rownames(mf) <- phy$tip.label
    mf[.order, ] <- tmp[1:nrow(tmp), ]
  }
  X <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  if (sum(!(y %in% c(0, 1)))) {
    stop("PGLMM.binary requires a binary response (dependent variable).")
  }
  if (var(y) == 0) {
    stop("The response (dependent variable) is always 0 or always 1.")
  }
  p <- ncol(X)
  Vphy <- ape::vcv(phy)
  Vphy <- Vphy/max(Vphy)
  Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)
  
  if (!is.null(B.init) & length(B.init) != p) {
    warning("B.init not correct length, so computed B.init using glm()")
  }
  if (is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) {
    B.init <- t(matrix(glm(formula = formula, data = data, family = "binomial")$coefficients, ncol = p))
  }
  B <- B.init
  s2 <- s2.init
  b <- matrix(0, nrow = n)
  beta <- rbind(B, b)
  mu <- exp(X %*% B)/(1 + exp(X %*% B))
  XX <- cbind(X, diag(1, nrow = n, ncol = n))
  C <- s2 * Vphy
  est.s2 <- s2
  est.B <- B
  oldest.s2 <- 10^6
  oldest.B <- matrix(10^6, nrow = length(est.B))
  iteration <- 0
  exitflag <- 0
  rcondflag <- 0
  while (((crossprod(est.s2 - oldest.s2) > tol.pql^2) | (crossprod(est.B - oldest.B)/length(B) > tol.pql^2)) & (iteration <= maxit.pql)) {
    iteration <- iteration + 1
    oldest.s2 <- est.s2
    oldest.B <- est.B
    est.B.m <- B
    oldest.B.m <- matrix(10^6, nrow = length(est.B))
    iteration.m <- 0
    while ((crossprod(est.B.m - oldest.B.m)/length(B) > tol.pql^2) & (iteration.m <= maxit.pql)) {
      iteration.m <- iteration.m + 1
      oldest.B.m <- est.B.m
      invW <- diag(as.vector((mu * (1 - mu))^-1))
      V <- invW + C
      if (sum(is.infinite(V)) > 0 | rcond(V) < 10^-10) {
        rcondflag <- rcondflag + 1
        B <- 0 * B.init + 0.001
        b <- matrix(0, nrow = n)
        beta <- rbind(B, b)
        mu <- exp(X %*% B)/(1 + exp(X %*% B))
        oldest.B.m <- matrix(10^6, nrow = length(est.B))
        invW <- diag(as.vector((mu * (1 - mu))^-1))
        V <- invW + C
      }
      invV <- solve(V)
      Z <- X %*% B + b + (y - mu)/(mu * (1 - mu))
      denom <- t(X) %*% invV %*% X
      num <- t(X) %*% invV %*% Z
      B <- as.matrix(solve(denom, num))
      b <- C %*% invV %*% (Z - X %*% B)
      beta <- rbind(B, b)
      mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))
      est.B.m <- B
    }
    H <- Z - X %*% B
    opt <- optim(fn = pglmm.reml, par = s2, tinvW = invW, tH = H, tVphy = Vphy, tX = X, method = "BFGS", control = list(factr = 1e+12, maxit = maxit.reml))
    s2 <- abs(opt$par)
    C <- s2 * Vphy
    est.s2 <- s2
    est.B <- B
  }
  convergeflag <- "converged"
  if (iteration >= maxit.pql | rcondflag >= 3) {
    convergeflag <- "Did not converge; try increasing maxit.pql or starting with B.init values of .001"
  }
  converge.test.s2 <- (crossprod(est.s2 - oldest.s2))^0.5
  converge.test.B <- (crossprod(est.B - oldest.B))^0.5/length(est.B)
  invW <- diag(as.vector((mu * (1 - mu))^-1))
  V <- invW + C
  invV <- solve(V)
  Z <- X %*% B + b + (y - mu)/(mu * (1 - mu))
  denom <- t(X) %*% invV %*% X
  num <- t(X) %*% invV %*% Z
  B <- solve(denom, num)
  b <- C %*% invV %*% (Z - X %*% B)
  beta <- rbind(B, b)
  mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))
  H <- Z - X %*% B
  B.cov <- solve(t(X) %*% invV %*% X)
  B.se <- as.matrix(diag(B.cov))^0.5
  B.zscore <- B/B.se
  B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)
  LL <- opt$value
  lnlike.cond.reml <- -0.5 * (n - p) * log(2 * pi) + 0.5 * determinant(t(X) %*% X)$modulus[1] - 0.5 * LL
  LL0 <- pglmm.reml(par = 0, tinvW = invW, tH = H, tVphy = Vphy, tX = X)
  lnlike.cond.reml0 <- -0.5 * (n - p) * log(2 * pi) + 0.5 * determinant(t(X) %*% X)$modulus[1] - 0.5 * LL0
  P.H0.s2 <- pchisq(2 * (lnlike.cond.reml - lnlike.cond.reml0), df = 1, lower.tail = F)/2
  
  results <- list(formula = formula, B = B, B.se = B.se, B.cov = B.cov, B.zscore = B.zscore, B.pvalue = B.pvalue, s2 = s2, P.H0.s2 = P.H0.s2, mu = mu, b = b, B.init = B.init, 
    X = X, y = y, phy = phy, data = data, H = H, VCV = Vphy, V = V, convergeflag = convergeflag, iteration = iteration, converge.test.s2 = converge.test.s2, converge.test.B = converge.test.B, 
    rcondflag = rcondflag)
  class(results) <- "binaryPGLMM"
  
  results
}
