#' Calculate R2.ce
#'
#' Calculate R2.ce for LMM, GLMM, PGLM, and PGLMMs.
#' @param mod a regression model with the following class: 'lmerMod', 'glmerMod', 'phylolm', and 'binaryPGLMM'
#' @param mod.r reduced model, if not provided, will use corresponding models with intercept as the only predictor
#' @param phy a phylogeny with tip labels and branch length
#' @param cpp whether to use Rcpp version of loops when the mod is class of 'phylolm', default is TRUE
#' @return R2.ce
#' @export
#'
R2.ce <- function(mod = NULL, mod.r = NULL, phy = NULL, cpp = TRUE) {
    
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
        return(R2.ce.lmerMod(mod, mod.r))
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
        return(R2.ce.glmerMod(mod, mod.r)[1])
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
        if (cpp) {
            return(R2.ce.phylolm(mod, mod.r, phy, cpp = TRUE))
        } else {
            return(R2.ce.phylolm(mod, mod.r, phy, cpp = FALSE))
        }
    }
    
    if (class(mod)[1] == "binaryPGLMM") {
        if (!is.object(mod.r)) {
            Y <- mod$y
            mod.r <- glm(Y ~ 1, family = "binomial")
        }
        if (!is.element(class(mod.r)[1], c("binaryPGLMM", "glm"))) {
            stop("mod.r must be class binaryPGLMM or glm.")
        }
        return(R2.ce.binaryPGLMM(mod, mod.r)[1])
    }
}

R2.ce.lmerMod <- function(mod = NULL, mod.r = NULL) {
    Y <- model.frame(mod)[, 1]
    SSE.ce <- var(Y - stats::fitted(mod))
    SSE.ce.r <- var(Y - stats::fitted(mod.r))
    return(1 - SSE.ce/SSE.ce.r)
}

R2.ce.glmerMod <- function(mod = NULL, mod.r = NULL) {
    Y <- model.frame(mod)[, 1]
    SSE.ce <- var(Y - stats::fitted(mod))
    SSE.ce.r <- var(Y - stats::fitted(mod.r))
    return(1 - SSE.ce/SSE.ce.r)
}

R2.ce.phylolm <- function(mod = NULL, mod.r = NULL, phy = NULL, cpp = TRUE) {
    Y <- mod$y
    X <- mod$X
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    optpar <- round(mod$optpar, digits = 4)
    if (mod$model == "lambda") {
        phy.f <- phylolm::transf.branch.lengths(phy, parameters = list(lambda = optpar), model = mod$model)$tree
    }
    if (mod$model %in% c("OUrandomRoot", "OUfixedRoot")) {
        phy.f <- phylolm::transf.branch.lengths(phy, parameters = list(alpha = optpar), model = mod$model)$tree
    }
    
    V <- ape::vcv(phy.f)
    R <- Y - stats::fitted(mod)
    
    if (cpp) {
        Rhat = loop_cpp(as.matrix(R), as.matrix(V))
    } else {
        Rhat <- matrix(0, nrow = n, ncol = 1)
        for (j in 1:n) {
            r <- R[-j]
            VV <- V[-j, -j]
            iVV <- solve(VV)
            v <- V[j, -j]
            Rhat[j] <- v %*% iVV %*% r
        }
    }
    
    Yhat <- as.numeric(stats::fitted(mod) + Rhat)
    SSE.ce <- var(Y - Yhat)
    
    # reduced model
    if (class(mod.r) == "phylolm") {
        X.r <- mod.r$X
        optpar.r <- round(mod.r$optpar, digits = 4)
        if (mod$model == "lambda") {
            phy.r <- phylolm::transf.branch.lengths(phy, parameters = list(lambda = optpar.r), model = mod.r$model)$tree
        }
        if (mod$model %in% c("OUrandomRoot", "OUfixedRoot")) {
            phy.r <- phylolm::transf.branch.lengths(phy, parameters = list(alpha = optpar.r), model = mod.r$model)$tree
        }
        V.r <- ape::vcv(phy.r)
        R.r <- Y - stats::fitted(mod.r)
        if (cpp) {
            Rhat.r = loop_cpp(as.matrix(R.r), as.matrix(V.r))
        } else {
            Rhat.r <- matrix(0, nrow = n, ncol = 1)
            for (j in 1:n) {
                r.r <- R.r[-j]
                VV.r <- V.r[-j, -j]
                iVV.r <- solve(VV.r)
                v.r <- V.r[j, -j]
                Rhat.r[j] <- v.r %*% iVV.r %*% r.r
            }
        }
        Yhat.r <- as.numeric(stats::fitted(mod.r) + Rhat.r)
    }
    
    if (class(mod.r) == "lm") {
        Yhat.r <- stats::fitted(mod.r)
    }
    
    SSE.ce.r <- var(Y - Yhat.r)
    
    return(1 - SSE.ce/SSE.ce.r)
}

R2.ce.binaryPGLMM <- function(mod = NULL, mod.r = NULL) {
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
    
    return(1 - SSE.ce/SSE.ce.r)
}
