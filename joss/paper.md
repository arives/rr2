---
title: '`rr2`: An R package to calculate $R^2$s for regression models'
tags:
  - R
  - $R^2$
  - GLMM
  - phylogenetic regression
  - models with correlated data
authors:
  - name: Anthony R. Ives
    orcid: 0000-0001-9375-9523
    affiliation: 1
  - name: Daijiang Li
    orcid: 0000-0002-0925-3421
    affiliation: 2
affiliations:
 - name: Department of Integrative Biology, UW-Madison, Madison, WI 53706
   index: 1
 - name: Department of Wildlife Ecology and Conservation, University of Florida, Gainesville, FL 32611
   index: 2
date: 6 September 2018
bibliography: paper.bib
---

# Summary

Reporting the variance explained by a model (an $R^2$) is common for many simple statistical tests. However, conceptual challenges exist in defining $R^2$ for models that include correlated data. @ives2018r2s proposed three $R^2$s ($R^2_{lik}$, $R^2_{resid}$, and $R^2_{pred}$) for a variety of regression models that include correlation among data such as linear mixed models (LMMs), generalized linear mixed models (GLMMs), and phylogenetic regressions [PGLMMs, @ives2011generalized; @ives2014phylogenetic]. These three $R^2$s can also be used as partial $R^2$s to compare the contributions of predictor variables (fixed effects) and/or correlation structures (random effects) to the fit of the models.

The `rr2` package provides R functions to implement the $R^2$s propsed by @ives2018r2s. The main function is `R2()`, which calculates all three $R^2$s by default, with arguments available to select which $R^2$(s) to calculate by users. Alternatively, individual $R^2$s can be calculated with corresponding functions (`R2_lik()`, `R2_resid()`, and `R2_pred()`). Supported models include linear models (`lm`), generalized linear models (`glm`), linear mixed models (`lmerMod`), generalized linear mixed models (`glmerMod`), phylogenetic linear models (`phylolm`), and phylogenetic generalized linear mixed models (`binaryPGLMM`, `phyloglm`, and `communityPGLMM`). 

The R package `rr2` is available on [Github](https://github.com/arives/rr2), where issues can be opened.

# Acknowledgments

This work was funded by NSF grants NSF/NASA-DEB-Dimensions 1240804 and DEB-LTREB-1052160.

# References
