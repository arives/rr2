#' Calculate R2.lik, R2.resid, and R2.pred
#'
#' Calculate all three R2s for LMM, GLMM, PGLM, and PGLMMs.
#' @param mod a regression model with the following class: 'lmerMod', 'glmerMod', 'phylolm', and 'binaryPGLMM'
#' @param mod.r reduced model, if not provided, will use corresponding models with intercept as the only predictor
#' @param phy the phylogeny for phylogenetic models, which is not required to be specified for R2.lik.
#' @param lik whether to calculate R2.lik, default is TRUE
#' @param resid whether to calculate R2.resid, default is TRUE
#' @param pred whether to calculate R2.pred, default is TRUE
#' @return all three R2s
#' @export
#'
R2 = function(mod = NULL, mod.r = NULL, phy = NULL, lik = TRUE, resid = TRUE, pred = TRUE){
  
  # phyloglm only have R2.lik method
  if(class(mod) == "phyloglm"){
    resid = FALSE; pred = FALSE
    message("models with class phyloglm only have R2.lik method")
  }
  
  # binaryPGLMM does not have R2.lik method
  if(class(mod) == "binaryPGLMM"){
    lik = FALSE
    message("models with class binaryPGLMM do not have R2.lik method")
  }
  
  # phylolm requires phy object
  if(class(mod) == "phylolm" & is.null(phy)) stop("phy object is required for models with class phylolm")
  
  out = data.frame(R2s = c("R2_lik", "R2_resid", "R2_pred"), value = NA, stringsAsFactors = FALSE)
  
  if(is.null(phy)){
    if(lik) out$value[1] = R2.lik(mod, mod.r)
    if(resid) out$value[2] = R2.resid(mod, mod.r)
    if(pred) out$value[3] = R2.pred(mod, mod.r)
  } else {
    if(lik) out$value[1] = R2.lik(mod, mod.r)
    if(resid) out$value[2] = R2.resid(mod, mod.r, phy)
    if(pred) out$value[3] = R2.pred(mod, mod.r, phy)
  }
  
  out = na.omit(out) # remove R2s not calculated
  row.names(out) = NULL # reset row names
  if(nrow(out) == 0) warning("at least to calculate one R2")
  
  return(out)
} 
