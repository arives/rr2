#' Calculate R2.lik, R2.resid, and R2.pred
#'
#' This is a wrapper for calculating all three R2s -- `R2.lik`, `R2.resid`, and 
#' `R2.pred` -- for LMM, GLMM, PGLM, and PGLMMs. Note that the individual functions
#' `R2.lik`, `R2.resid`, and `R2.pred` can be called. This is preferrable if you are only 
#' interested in one R2; for example, for 'phylolm' called from `R2` you need to specify 
#' the phy (phylo object for the phylogeny), while `R2.lik` does not require this.
#'   
#' @param mod A regression model with the following class: 'lmerMod', 'glmerMod', 'phylolm', 'binaryPGLMM', and 'communityPGLMM'.
#' @param mod.r A reduced model, if not provided, will use corresponding models with intercept as the only predictor.
#' @param phy The phylogeny for phylogenetic models, which is not required to be specified for R2.lik.
#' @param sigma2_d Distribution-specific variance for logistic regressions (GLM, GLMM, and PGLMM). 
#'    Available options are 'corrected' and 'NS' (stands for Nakagawa and Schielzeth 2013).
#' @param lik Whether to calculate R2.lik, default is TRUE.
#' @param resid Whether to calculate R2.resid, default is TRUE.
#' @param pred Whether to calculate R2.pred, default is TRUE.
#' @return A vector, with all three R2s by default.
#' @export
#'
R2 <- function(mod = NULL, mod.r = NULL, phy = NULL, sigma2_d = c('corrected', 'NS'), 
               lik = TRUE, resid = TRUE, pred = TRUE) {
  
  if(all(!lik, !resid, !pred)) 
    stop("Specify at least one of 'lik', 'resid', or 'pred' for R2")
  
  sigma2_d <- match.arg(sigma2_d)
  
  # phyloglm only have R2.lik method
  if (any(class(mod) %in% "phyloglm")) {
    resid <- FALSE
    pred <- FALSE
    message("Models of class phyloglm only have R2.lik method")
  }
  
  # gaussian communityPGLMM only have R2.lik method
  if (any(class(mod) %in% "communityPGLMM")) {
    if(mod$family == "gaussian"){
      resid <- FALSE
      message("Models of class communityPGLMM (gaussian) do not have R2.resid method")
    }
  }
  
  # binaryPGLMM does not have R2.lik method
  if (any(class(mod) %in% "binaryPGLMM") & lik == TRUE) {
    lik <- FALSE
    message("Models of class binaryPGLMM do not have R2.lik method")
  }
  
  # binary communityPGLMM only have R2.pred method at this moment
  if (any(class(mod) %in% "communityPGLMM") & (lik == TRUE | resid == TRUE)) {
    if(mod$family == "binomial"){
      resid <- FALSE
      lik <- FALSE
      message("Models of class communityPGLMM (binomial) only have R2.pred method")
    }
  }
  
  # phylolm requires phy object except for R2.lik
  if (any(class(mod) %in% "phylolm") & is.null(phy)) 
    stop("Phy object is required for models with class phylolm")
  
  out <- rep(NA, 3)
  names(out) <- c("R2_lik", "R2_resid", "R2_pred")
  if (lik) 
    out[1] <- R2.lik(mod, mod.r)
  if (resid) 
    out[2] <- R2.resid(mod, mod.r, phy, sigma2_d)
  if (pred) 
    out[3] <- R2.pred(mod, mod.r, phy)
  
  out <- out[!is.na(out)] # remove R2s not calculated
  
  return(out)
}
