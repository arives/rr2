#' Calculate R2.lr, R2.ls, and R2.ce
#'
#' Calculate all three R2s for LMM, GLMM, PGLM, and PGLMMs.
#' @param mod a regression model with the following class: 'lmerMod', 'glmerMod', 'phylolm', and 'binaryPGLMM'
#' @param mod.r reduced model, if not provided, will use corresponding models with intercept as the only predictor
#' @param phy the phylogeny for phylogenetic models, which is not required to be specified for R2.lr.
#' @param lr whether to calculate R2.lr, default is TRUE
#' @param ls whether to calculate R2.ls, default is TRUE
#' @param ce whether to calculate R2.ce, default is TRUE
#' @return all three R2s
#' @export
#'
R2 = function(mod = NULL, mod.r = NULL, phy = NULL, lr = TRUE, ls = TRUE, ce = TRUE){
  
  out = data.frame(R2s = c("R2_lr", "R2_ls", "R2_ce"), value = NA)
  
  if(is.null(phy)){
    if(lr) out$value[1] = R2.lr(mod, mod.r)
    if(ls) out$value[2] = R2.ls(mod, mod.r)
    if(ce) out$value[3] = R2.ce(mod, mod.r)
  } else {
    if(lr) out$value[1] = R2.lr(mod, mod.r)
    if(ls) out$value[2] = R2.ls(mod, mod.r, phy)
    if(ce) out$value[3] = R2.ce(mod, mod.r, phy)
  }
  
  out = na.omit(out) # remove R2s not calculated
  row.names(out) = NULL # reset row names
  if(nrow(out) == 0) warning("at least to calculate one R2")
  
  return(out)
} 
