# R2_source_code

require(lme4)
require(phylolm)
# require(nlme)
# require(ape)

########################################################
########################################################
# utility functions not used in primary functions
########################################################
########################################################

inv.logit <- function(x){
	1/(1 + exp(-x))
}

partialR2 <- function(mod, mod.r){
    anova.full <- anova(mod)
    anova.reduced <- anova(mod.r)

    sse.full <- tail(anova.full$"Sum Sq", 1)
    sse.reduced <- tail(anova.reduced$"Sum Sq", 1)

    pR2 <- (sse.reduced - sse.full) / sse.reduced
    return(pR2)
}    

partialR2adj <- function(mod, df.f, mod.r, df.r){
    anova.f <- anova(mod)
    anova.r <- anova(mod.r)

    sse.f <- tail(anova.f$"Sum Sq", 1)
    sse.r <- tail(anova.r$"Sum Sq", 1)

    R2 <- 1 - sse.f / sse.r
    R2.adj <- 1 - (sse.f/df.f) / (sse.r/df.r)
    return(list(R2 = R2, R2.adj = R2.adj))
}  
  
########################################################
########################################################
# R2.ls
########################################################
########################################################

R2.ls <- function(mod = NA, mod.r = NA, phy = NA) {
	
	if(!is.element(class(mod)[1], c("lmerMod","glmerMod","phylolm","binaryPGLMM"))) {
		stop("mod must be class one of classes lmerMod, glmerMod, phylolm, binaryPGLMM.")
	}
	if(class(mod)[1] == "lmerMod"){
		if(!is.object(mod.r)){
			Y <- model.frame(mod)[,1]
			mod.r <- lm(Y ~ 1)
		}
		if(!is.element(class(mod)[1], c("lmerMod","lm"))){
			stop("mod.r must be class lmerMod or lm.")
		}
		return(R2.ls.lmerMod(mod, mod.r))
	}
	if(class(mod)[1] == "glmerMod"){
		if(family(mod)[[1]] != "binomial") {
			stop("Sorry, but only binomial (binary) glmerMod models are allowed.")
		}
		if(!is.object(mod.r)){
			Y <- model.frame(mod)[,1]
			mod.r <- glm(Y ~ 1, family=family(mod)[[1]])
		}
		if(!is.element(class(mod.r)[1], c("glmerMod","glm"))){
			stop("mod.r must be class glmerMod or glm.")
		}
		return(R2.ls.glmerMod(mod, mod.r)[1])
	}
	if(class(mod)[1] == "phylolm"){
		if(!is.object(phy)){
			stop("For phylolm you must provide the phylo object")
		}
		if(!is.object(mod.r)){
			Y <- mod$y
			mod.r <- lm(Y ~ 1)
		}
		if(!is.element(class(mod.r)[1], c("phylolm","lm"))){
			stop("mod.r must be class phylolm or lm.")
		}
		return(R2.ls.phylolm(mod, mod.r, phy))
	}
	if(class(mod)[1] == "binaryPGLMM"){
		if(!is.object(mod.r)){
			Y <- mod$y
			mod.r <- glm(Y ~ 1, family = "binomial")
		}
		if(!is.element(class(mod.r)[1], c("binaryPGLMM","glm"))){
			stop("mod.r must be class binaryPGLMM or glm.")
		}
		return(R2.ls.binaryPGLMM(mod, mod.r)[1])
	}
}

R2.ls.lmerMod <- function(mod = NA, mod.r = NA) {
	
	Y <- model.frame(mod)[,1]
	X <- model.matrix(mod)
	n <- dim(X)[1]
	p <- dim(X)[2]
	
	X.r <- model.matrix(mod.r)
	p.r <- dim(X.r)[2]

	D <- attr(mod,"pp")$LamtUt
	V <- t(D) %*% D + diag(n)
	
	if(class(mod.r) == "lmerMod"){
		D.r <- attr(mod.r,"pp")$LamtUt
		V.r <- t(D.r) %*% D.r + diag(n)
	}
	if(class(mod.r) == "lm"){
		V.r <- diag(n)
	}

	vcov <- as.data.frame(VarCorr(mod))$vcov
	sigma2 <- vcov[length(vcov)]
	if(class(mod.r) == "lmerMod"){
		vcov.r <- as.data.frame(VarCorr(mod.r))$vcov
		sigma2.r <- vcov.r[length(vcov.r)]
	}
	if(class(mod.r) == "lm"){
		sigma2.r <- (n-p.r)/n * sigma(mod.r)^2
	}
	R2.ls <- 1 - sigma2/sigma2.r
	return(R2.ls)
}

R2.ls.glmerMod <- function(mod = NA, mod.r = NA) {
	
	Y <- model.frame(mod)[,1]
	X <- model.matrix(mod)
	n <- dim(X)[1]
	p <- dim(X)[2]
	
	X.r <- model.matrix(mod.r)
	p.r <- dim(X.r)[2]
	
	# full model
	AVAg <- t(attr(mod,"pp")$"LamtUt") %*% attr(mod,"pp")$"LamtUt"
	Ag <- diag(1/attr(mod,"pp")$"Xwts")
	
	mu <- family(mod)$linkinv(X %*% fixef(mod))

	g <- mu*(1-mu)
	g <- g/prod(g)^(1/n)
	V <- Ag %*% AVAg %*% Ag + diag(n)
	iV <- solve(V)

	iA <- diag(as.numeric(g^(-1/2)))
	iAVA <- iA %*% iV %*% iA	
	SSE.ls <- t(Y - mu) %*% iAVA %*% (Y - mu)

	# reduced model
	if(class(mod.r)[1] == "glmerMod"){
		AVAg.r <- t(attr(mod.r,"pp")$"LamtUt") %*% attr(mod.r,"pp")$"LamtUt"
		Ag.r <- diag(1/attr(mod.r,"pp")$"Xwts")
		
		mu.r <- family(mod.r)$linkinv(X.r %*% fixef(mod.r))
		g.r <- mu.r*(1-mu.r)
		g.r <- g.r/prod(g.r)^(1/n)
		V.r <- Ag.r %*% AVAg.r %*% Ag.r + diag(n)
		iV.r <- solve(V.r)
	
		iA.r <- diag(as.numeric(g.r^(-1/2)))
		iAVA.r <- iA.r %*% iV.r %*% iA.r	
		SSE.ls.r <- t(Y - mu.r) %*% iAVA.r %*% (Y - mu.r)
	}
	if(class(mod.r)[1] == "glm"){
		mu.r <- mod.r$fitted.values
		g.r <- mu.r*(1-mu.r)
		g.r <- g.r/prod(g.r)^(1/n)
		V.r <- diag(g.r)
		iV.r <- diag(1/g.r)
		SSE.ls.r <- t(Y - mu.r) %*% iV.r %*% (Y - mu.r)
	}
	R2.ls <- 1 - SSE.ls/SSE.ls.r
	return(R2.ls)
}

R2.ls.phylolm <- function(mod = NA, mod.r = NA, phy = NA) {
	
	Y <- mod$y
	X <- mod$X
	n <- dim(X)[1]
	p <- dim(X)[2]
	
	optpar <- round(mod$optpar, digits=4)
	if(mod$model == "lambda") phy.f <- transf.branch.lengths(phy, parameters = list("lambda"=optpar), model = mod$model)$tree
	if(mod$model == "OUrandomRoot") phy.f <- transf.branch.lengths(phy, parameters = list("alpha"=optpar), model = mod$model)$tree
	if(mod$model == "OUfixedRoot") phy.f <- transf.branch.lengths(phy, parameters = list("alpha"=optpar), model = mod$model)$tree
	V <- vcv(phy.f)
	scal <- sum(phy.f$edge.length)/n
	V <- V/scal
	
	sigma2 <- mod$sigma2

	if(class(mod.r) == "phylolm"){
		X.r <- mod.r$X
		p.r <- dim(X.r)[2]

		optpar.r <- round(mod.r$optpar, digits=4)
		if(mod$model == "lambda") phy.r <- transf.branch.lengths(phy, parameters = list("lambda"=optpar.r), model = mod.r$model)$tree
		if(mod$model == "OUrandomRoot") phy.r <- transf.branch.lengths(phy, parameters = list("alpha"=optpar.r), model = mod.r$model)$tree
		if(mod$model == "OUfixedRoot") phy.r <- transf.branch.lengths(phy, parameters = list("alpha"=optpar.r), model = mod.r$model)$tree
		V.r <- vcv(phy.r)
		scal.r <- sum(phy.r$edge.length)/n
		V.r <- V.r/scal.r
		sigma2.r <- mod.r$sigma2
	}
	if(class(mod.r) == "lm"){
		X.r <- model.matrix(mod.r)
		p.r <- dim(X.r)[2]
		V.r <- diag(n)
		scal.r <- 1
		sigma2.r <- (n-p.r)/n * sigma(mod.r)^2
	}

	R2.ls <- 1 - (scal*sigma2)/(scal.r*sigma2.r)
	return(R2.ls)
}

R2.ls.binaryPGLMM <- function(mod = NA, mod.r = NA) {
	
	# utility function
	inv.logit <- function(x){
		1/(1 + exp(-x))
	}

	mu <- mod$mu
	n <- length(mu)
	Yhat <- mu
	R <- mod$y - mu
	
	if(class(mod.r)[1] == "binaryPGLMM") {
		mu.r <- mod.r$mu
		Yhat.r <- mu.r
		R.r <- mod.r$y - mu.r
	}else{
		mu.r <- mod.r$fitted.values
		Yhat.r <- mu.r
		R.r <- mod.r$y - mu.r
	}
	
	phy <- mod$phy
	
	g <- mu*(1-mu)
	g <- g/prod(g)^(1/n)
	V <- mod$s2 * vcv(phy) + diag(n)
	phy.temp <- vcv2phylo(V)
	V <- V/(sum(phy.temp$edge.length)/n)	
	iV <- solve(V)
	iA <- diag(as.numeric(g^(-1/2)))
	
	iAVA <- iA %*% iV %*% iA	
	SSE.ls <- t(R) %*% iAVA %*% R

	# reduced model
	if(class(mod.r)[1] == "binaryPGLMM") {
		g.r <- mu.r*(1-mu.r)
		g.r <- g.r/prod(g.r)^(1/n)
		V.r <- mod.r$s2 * vcv(phy) + diag(n)
		phy.temp <- vcv2phylo(V.r)
		V.r <- V.r/(sum(phy.temp$edge.length)/n)
		iV.r <- solve(V.r)
		iA.r <- diag(as.numeric(g.r^(-1/2)))
		
		iAVA.r <- iA.r %*% iV.r %*% iA.r	
		SSE.ls.r <- t(R.r) %*% iAVA.r %*% R.r
	}
	if(class(mod.r)[1] == "glm"){
		mu.r <- mod.r$fitted.values
		g.r <- mu.r*(1-mu.r)
		g.r <- g.r/prod(g.r)^(1/n)
		V.r <- diag(g.r)
		iV.r <- diag(1/g.r)
		SSE.ls.r <- t(R.r) %*% iV.r %*% R.r
	}
	R2.ls <- 1 - SSE.ls/SSE.ls.r
	return(R2.ls)
}


########################################################
########################################################
# R2.ce
########################################################
########################################################

R2.ce <- function(mod = NA, mod.r = NA, phy = NA) {
	
	if(!is.element(class(mod)[1], c("lmerMod","glmerMod","phylolm","binaryPGLMM"))) {
		stop("mod must be class one of classes lmerMod, glmerMod, phylolm, binaryPGLMM.")
	}
	if(class(mod)[1] == "lmerMod"){
		if(!is.object(mod.r)){
			Y <- model.frame(mod)[,1]
			mod.r <- lm(Y ~ 1)
		}
		if(!is.element(class(mod)[1], c("lmerMod","lm"))){
			stop("mod.r must be class lmerMod or lm.")
		}
		return(R2.ce.lmerMod(mod, mod.r))
	}
	if(class(mod)[1] == "glmerMod"){
		if(family(mod)[[1]] != "binomial") {
			stop("Sorry, but only binomial (binary) glmerMod models are allowed.")
		}
		if(!is.object(mod.r)){
			Y <- model.frame(mod)[,1]
			mod.r <- glm(Y ~ 1, family=family(mod)[[1]])
		}
		if(!is.element(class(mod.r)[1], c("glmerMod","glm"))){
			stop("mod.r must be class glmerMod or glm.")
		}
		return(R2.ce.glmerMod(mod, mod.r)[1])
	}
	if(class(mod)[1] == "phylolm"){
		if(!is.object(phy)){
			stop("For phylolm you must provide the phylo object")
		}
		if(!is.object(mod.r)){
			Y <- mod$y
			mod.r <- lm(Y ~ 1)
		}
		if(!is.element(class(mod.r)[1], c("phylolm","lm"))){
			stop("mod.r must be class phylolm or lm.")
		}
		return(R2.ce.phylolm(mod, mod.r, phy))
	}
	if(class(mod)[1] == "binaryPGLMM"){
		if(!is.object(mod.r)){
			Y <- mod$y
			mod.r <- glm(Y ~ 1, family = "binomial")
		}
		if(!is.element(class(mod.r)[1], c("binaryPGLMM","glm"))){
			stop("mod.r must be class binaryPGLMM or glm.")
		}
		return(R2.ce.binaryPGLMM(mod, mod.r)[1])
	}
}

R2.ce.lmerMod <- function(mod = NA, mod.r = NA) {
	
	Y <- model.frame(mod)[,1]
	X <- model.matrix(mod)
	n <- dim(X)[1]
	p <- dim(X)[2]
	
	X.r <- model.matrix(mod.r)
	p.r <- dim(X.r)[2]

	D <- attr(mod,"pp")$LamtUt
	V <- t(D) %*% D + diag(n)
	
	if(class(mod.r) == "lmerMod"){
		D.r <- attr(mod.r,"pp")$LamtUt
		V.r <- t(D.r) %*% D.r + diag(n)
	}
	if(class(mod.r) == "lm"){
		V.r <- diag(n)
	}

	iV <- solve(V)
		
	bhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% Y)
	fitted.values <- X %*% bhat
	R <- Y - X %*% bhat

	Rhat <- matrix(0, nrow=n, ncol=1)
	for(j in 1:n){
		
		r <- R[-j]
		VV <- V[-j,-j]
		iVV <- solve(VV)

		# version using global mean: This works best
		bbhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% R)[1]		
		v <- V[j, -j]		
		Rhat[j] <- as.numeric(bbhat + v %*% iVV %*% (r - bbhat))
	}
	Yhat <- as.numeric(fitted.values + Rhat)
	SSE.ce <- var(Y-Yhat)
	
	# reduced model
	iV.r <- solve(V.r)
		
	bhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV.r %*% Y)
	fitted.values <- X.r %*% bhat
	R <- Y - X.r %*% bhat

	Rhat <- matrix(0, nrow=n, ncol=1)
	for(j in 1:n){
		
		r <- R[-j]
		VV.r <- V.r[-j,-j]
		iVV.r <- solve(VV.r)

		# version using global mean: This works best
		bbhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV %*% R)[1]		
		v <- V.r[j, -j]		
		Rhat[j] <- as.numeric(bbhat + v %*% iVV.r %*% (r - bbhat))
	}
	Yhat.r <- as.numeric(fitted.values + Rhat)
	SSE.ce.r <- var(Y-Yhat.r)

	R2.ce <- 1 - SSE.ce/SSE.ce.r
	return(R2.ce)
}

R2.ce.glmerMod <- function(mod = NA, mod.r = NA) {
	
	Y <- model.frame(mod)[,1]
	X <- model.matrix(mod)

	n <- dim(X)[1]
	p <- dim(X)[2]
	
	X.r <- model.matrix(mod.r)
	p.r <- dim(X.r)[2]
	
	AVAg <- t(attr(mod,"pp")$"LamtUt") %*% attr(mod,"pp")$"LamtUt"
	Ag <- diag(1/attr(mod,"pp")$"Xwts")
	
	mu <- family(mod)$linkinv(X %*% fixef(mod))
	mu[mu<10^-6] <- 10^-6
	mu[mu>(1-10^-6)] <- (1-10^-6)

	g <- mu*(1-mu)
#	g <- g/prod(g)^(1/n)
	V <- Ag %*% AVAg %*% Ag + diag(n)
	iV <- solve(V)

	iA <- diag(as.numeric(g^(-1/2)))
	iAVA <- iA %*% iV %*% iA	

	A <- diag(as.numeric(g^(1/2)))
	AVA <- A %*% V %*% A	

	Yhat <- mu

	# reduced model
	if(class(mod.r)[1] == "glmerMod"){
		AVAg.r <- t(attr(mod.r,"pp")$"LamtUt") %*% attr(mod.r,"pp")$"LamtUt"
		Ag.r <- diag(1/attr(mod.r,"pp")$"Xwts")
		
		mu.r <- family(mod.r)$linkinv(X.r %*% fixef(mod.r))
		mu.r[mu.r<10^-6] <- 10^-6
		mu.r[mu.r>(1-10^-6)] <- (1-10^-6)

		g.r <- mu.r*(1-mu.r)
#		g.r <- g.r/prod(g.r)^(1/n)
		V.r <- Ag.r %*% AVAg.r %*% Ag.r + diag(n)
		iV.r <- solve(V.r)
	
		iA.r <- diag(as.numeric(g.r^(-1/2)))
		iAVA.r <- iA.r %*% iV.r %*% iA.r	

		A.r <- diag(as.numeric(g.r^(1/2)))
		AVA.r <- A.r %*% V.r %*% A.r	

		Yhat.r <- mu.r
	}
	if(class(mod.r)[1] == "glm"){
		mu.r <- mod.r$fitted.values
		mu.r[mu.r<10^-6] <- 10^-6
		mu.r[mu.r>(1-10^-6)] <- (1-10^-6)

		g.r <- mu.r*(1-mu.r)
#		g.r <- g.r/prod(g.r)^(1/n)
		AVA.r <- diag(g.r)
		iAVA.r <- diag(1/g.r)
		V.r <- diag(n)
		iV.r <- diag(n)

		Yhat.r <- mu.r
	}

	R <- (Y - Yhat)
	X <- matrix(1, ncol=1, nrow=n)
#	bbhat <- solve(t(X) %*% iAVA %*% X, t(X) %*% iAVA %*% R)[1]
	bbhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% R)[1]

	Rhat <- matrix(0, nrow=n, ncol=1)
	for(j in 1:n){
		
		r <- R[-j]
#		VV <- AVA[-j,-j]
		VV <- V[-j,-j]
		iVV <- solve(VV)
		
#		v <- AVA[j, -j]		
		v <- V[j, -j]		
		Rhat[j] <- as.numeric(bbhat + v %*% iVV %*% (r - bbhat))
	}
	Yhat.f <- Yhat + Rhat
	SSE.ce <- var(Y-Yhat.f)
	
	# reduced model	
	R.r <- (Y - Yhat.r)
	X <- matrix(1, ncol=1, nrow=n)
#	bbhat.r <- solve(t(X) %*% iAVA.r %*% X, t(X) %*% iAVA.r %*% R.r)[1]
	bbhat.r <- solve(t(X) %*% iV.r %*% X, t(X) %*% iV.r %*% R.r)[1]

	Rhat.r <- matrix(0, nrow=n, ncol=1)
	for(j in 1:n){
		
		r.r <- R.r[-j]
#		VV.r <- AVA.r[-j,-j]
		VV.r <- V.r[-j,-j]
		iVV.r <- solve(VV.r)
		
#		v.r <- AVA.r[j, -j]		
		v.r <- V.r[j, -j]		
		Rhat.r[j] <- as.numeric(bbhat.r + v.r %*% iVV.r %*% (r.r - bbhat.r))
	}
	Yhat.r <- Yhat.r + Rhat.r
	SSE.ce.r <- var(Y-Yhat.r)

	R2.ce <- 1 - SSE.ce/SSE.ce.r
	return(R2.ce)
}


R2.ce.phylolm <- function(mod = NA, mod.r = NA, phy = NA) {
	
	Y <- mod$y
	X <- mod$X
	n <- dim(X)[1]
	p <- dim(X)[2]
	
	optpar <- round(mod$optpar, digits=4)
	if(mod$model == "lambda") phy.f <- transf.branch.lengths(phy, parameters = list("lambda"=optpar), model = mod$model)$tree
	if(mod$model == "OUrandomRoot") phy.f <- transf.branch.lengths(phy, parameters = list("alpha"=optpar), model = mod$model)$tree
	if(mod$model == "OUfixedRoot") phy.f <- transf.branch.lengths(phy, parameters = list("alpha"=optpar), model = mod$model)$tree
	V <- vcv(phy.f)
	scal <- sum(phy.f$edge.length)/n
	V <- V/scal
	
	if(class(mod.r) == "phylolm"){
		X.r <- mod.r$X
		p.r <- dim(X.r)[2]

		optpar.r <- round(mod.r$optpar, digits=4)
		if(mod$model == "lambda") phy.r <- transf.branch.lengths(phy, parameters = list("lambda"=optpar.r), model = mod.r$model)$tree
		if(mod$model == "OUrandomRoot") phy.r <- transf.branch.lengths(phy, parameters = list("alpha"=optpar.r), model = mod.r$model)$tree
		if(mod$model == "OUfixedRoot") phy.r <- transf.branch.lengths(phy, parameters = list("alpha"=optpar.r), model = mod.r$model)$tree
		V.r <- vcv(phy.r)
		scal.r <- sum(phy.r$edge.length)/n
		V.r <- V.r/scal.r
	}
	if(class(mod.r) == "lm"){
		X.r <- model.matrix(mod.r)
		p.r <- dim(X.r)[2]
		V.r <- diag(n)
		scal.r <- 1
	}

	# full model
	iV <- solve(V)
		
	bhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% Y)
	fitted.values <- X %*% bhat
	R <- Y - X %*% bhat

	Rhat <- matrix(0, nrow=n, ncol=1)
	for(j in 1:n){
		
		r <- R[-j]
		VV <- V[-j,-j]
		iVV <- solve(VV)

		# version using global mean: This works best
		bbhat <- solve(t(X) %*% iV %*% X, t(X) %*% iV %*% R)[1]		
		v <- V[j, -j]		
		Rhat[j] <- bbhat + v %*% iVV %*% (r - bbhat)
	}
	Yhat <- as.numeric(fitted.values + Rhat)
	SSE.ce <- var(Y-Yhat)
	
	# reduced model
	iV.r <- solve(V.r)
		
	bhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV.r %*% Y)
	fitted.values <- X.r %*% bhat
	R <- Y - X.r %*% bhat

	Rhat <- matrix(0, nrow=n, ncol=1)
	for(j in 1:n){
		
		r <- R[-j]
		VV.r <- V.r[-j,-j]
		iVV.r <- solve(VV.r)

		# version using global mean: This works best
		bbhat <- solve(t(X.r) %*% iV.r %*% X.r, t(X.r) %*% iV %*% R)[1]		
		v <- V.r[j, -j]		
		Rhat[j] <- bbhat + v %*% iVV.r %*% (r - bbhat)
	}
	Yhat.r <- as.numeric(fitted.values + Rhat)
	SSE.ce.r <- var(Y-Yhat.r)

	R2.ce <- 1 - SSE.ce/SSE.ce.r
	return(R2.ce)
}

R2.ce.binaryPGLMM <- function(mod = NA, mod.r = NA) {
	
	# utility function
	inv.logit <- function(x){
		1/(1 + exp(-x))
	}

	mu <- mod$mu
	n <- length(mu)
	Yhat <- mu
	Y <- mod$y
	R <- Y - mu
	
	if(class(mod.r)[1] == "binaryPGLMM") {
		mu.r <- mod.r$mu
		Yhat.r <- mu.r
		Y <- mod.r$y
		R.r <- Y - mu.r
	}else{
		mu.r <- mod.r$fitted.values
		Yhat.r <- mu.r
		Y <- mod.r$y
		R.r <- Y - mu.r
	}
	
	phy <- mod$phy
	
	g <- mu*(1-mu)
	g <- g/prod(g)^(1/n)
	V <- mod$s2 * vcv(phy) + diag(n)
	phy.temp <- vcv2phylo(V)
	V <- V/(sum(phy.temp$edge.length)/n)	
	iV <- solve(V)

	iA <- diag(as.numeric(g^(-1/2)))
	iAVA <- iA %*% iV %*% iA	

	A <- diag(as.numeric(g^(1/2)))
	AVA <- A %*% V %*% A	

	Yhat <- mu

	# reduced model
	if(class(mod.r)[1] == "binaryPGLMM") {
		g.r <- mu.r*(1-mu.r)
		g.r <- g.r/prod(g.r)^(1/n)
		V.r <- mod.r$s2 * vcv(phy) + diag(n)
		phy.temp <- vcv2phylo(V.r)
		V.r <- V.r/(sum(phy.temp$edge.length)/n)
		iV.r <- solve(V.r)
	
		iA.r <- diag(as.numeric(g.r^(-1/2)))
		iAVA.r <- iA.r %*% iV.r %*% iA.r	
		
		A.r <- diag(as.numeric(g.r^(1/2)))		
		AVA.r <- A.r %*% V.r %*% A.r	
	}
	if(class(mod.r)[1] == "glm"){
		mu.r <- mod.r$fitted.values
		g.r <- mu.r*(1-mu.r)
		g.r <- g.r/prod(g.r)^(1/n)
		AVA.r <- diag(g.r)
		iAVA.r <- diag(1/g.r)	
	}
	
	R <- (Y - Yhat)
	X <- matrix(1, ncol=1, nrow=n)
	bbhat <- solve(t(X) %*% iAVA %*% X, t(X) %*% iAVA %*% R)[1]

	Rhat <- matrix(0, nrow=n, ncol=1)
	for(j in 1:n){
		
		r <- R[-j]
		VV <- AVA[-j,-j]
		iVV <- solve(VV)
		
		v <- AVA[j, -j]		
		Rhat[j] <- as.numeric(bbhat + v %*% iVV %*% (r - bbhat))
	}
	Yhat.f <- Yhat + Rhat
	SSE.ce <- var(Y-Yhat.f)
	
	# reduced model	
	R.r <- (Y - Yhat.r)
	X <- matrix(1, ncol=1, nrow=n)
	bbhat.r <- solve(t(X) %*% iAVA.r %*% X, t(X) %*% iAVA.r %*% R.r)[1]

	Rhat.r <- matrix(0, nrow=n, ncol=1)
	for(j in 1:n){
		
		r.r <- R.r[-j]
		VV.r <- AVA.r[-j,-j]
		iVV.r <- solve(VV.r)
		
		v.r <- AVA.r[j, -j]		
		Rhat.r[j] <- as.numeric(bbhat.r + v.r %*% iVV.r %*% (r.r - bbhat.r))
	}
	Yhat.r <- Yhat.r + Rhat.r
	SSE.ce.r <- var(Y-Yhat.r)

	R2.ce <- 1 - SSE.ce/SSE.ce.r
	return(R2.ce)
}

########################################################
########################################################
# R2.lr
########################################################
########################################################

R2.lr <- function(mod = NA, mod.r = NA) {
	
	if(!is.element(class(mod)[1], c("lmerMod","glmerMod","phylolm","phyloglm"))) {
		stop("mod must be class one of classes lmerMod, glmerMod, phylolm, phyloglm.")
	}
	if(class(mod)[1] == "lmerMod"){
		if(!is.object(mod.r)){
			Y <- model.frame(mod)[,1]
			mod.r <- lm(Y ~ 1)
		}
		if(!is.element(class(mod)[1], c("lmerMod","lm"))){
			stop("mod.r must be class lmerMod or lm.")
		}
		if(isREML(mod)){
			mod <- update(mod, REML=F)
			warning("mod updated with REML=F")
		}
		if(class(mod.r)[1] == "lmerMod" && isREML(mod.r)){
			mod.r <- update(mod.r, REML=F)
			warning("mod.r updated with REML=F")
		}
		return(R2.lr.lmerMod(mod, mod.r)[1])
	}
	if(class(mod)[1] == "glmerMod"){
		if(family(mod)[[1]] != "binomial") {
			stop("Sorry, but only binomial (binary) glmerMod models are allowed.")
		}
		if(!is.object(mod.r)){
			Y <- model.frame(mod)[,1]
			mod.r <- glm(Y ~ 1, family=family(mod)[[1]])
		}
		if(!is.element(class(mod.r)[1], c("glmerMod","glm"))){
			stop("mod.r must be class glmerMod or glm.")
		}
		return(R2.lr.glmerMod(mod, mod.r)[1])
	}
	if(class(mod)[1] == "phylolm"){
		if(!is.object(mod.r)){
			Y <- mod$y
			mod.r <- lm(Y ~ 1)
		}
		if(!is.element(class(mod.r)[1], c("phylolm","lm"))){
			stop("mod.r must be class phylolm or lm.")
		}
		return(R2.lr.phylolm(mod, mod.r))
	}
	if(class(mod)[1] == "phyloglm"){
		if(!is.object(mod.r)){
			Y <- mod$y
			mod.r <- glm(Y ~ 1, family="binomial")
		}
		if(!is.element(class(mod.r)[1], c("phyloglm","glm"))){
			stop("mod.r must be class phyloglm or glm.")
		}
		return(R2.lr.phyloglm(mod, mod.r)[1])
	}
}

R2.lr.lmerMod <- function(mod = NA, mod.r = NA) {
	
	X <- model.matrix(mod)
	n <- dim(X)[1]

	R2.lr <- 1-exp(-2/n * (logLik(mod) - logLik(mod.r)))
	return(R2.lr)
}

R2.lr.glmerMod <- function(mod = NA, mod.r = NA) {
	
	X <- model.matrix(mod)
	n <- dim(X)[1]

	R2.lr <- 1-exp(-2/n * (logLik(mod) - logLik(mod.r)))
	return(R2.lr)
}

R2.lr.phylolm <- function(mod = NA, mod.r = NA) {
	
	X <- mod$X
	n <- dim(X)[1]

	R2.lr <- 1-exp(-2/n * (logLik(mod)[[1]] - logLik(mod.r)[[1]]))
	return(R2.lr)
}

R2.lr.phyloglm <- function(mod = NA, mod.r = NA) {
	
	Y <- mod$y
	n <- dim(mod$X)[1]
	X <- mod$X[,-1]

	alpha.cutoff <- 40
	if(mod$alpha < alpha.cutoff) {
		LL <- mod$logLik
	}else{
		LL <- logLik(glm(Y ~ X, family="binomial"))
	}
	if(class(mod.r)[1] == "phyloglm"){
		if(mod.r$alpha < alpha.cutoff) {
			LL.r <- mod.r$logLik		
		}else{
			X.r <- mod.r$X
			LL.r <- logLik(glm(Y ~ 0 + X.r, family="binomial"))
		}
	}else{
		LL.r <- logLik(mod.r)
	}
	
	R2.lr <- (1 - exp(-2/n*(LL - LL.r)))/(1 - exp(2/n*LL.r))
	return(R2.lr)
}


########################################################
########################################################
# version of binaryPGLMM that returns values of Y
########################################################
########################################################

binaryPGLMM <- function (formula, data = list(), phy, s2.init = 0.1, B.init = NULL, tol.pql = 10^-6, maxit.pql = 200, maxit.reml = 100) 
{
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
                LL <- logdetV + t(tH) %*% invV %*% tH + determinant(t(tX) %*% 
                  invV %*% tX)$modulus[1]
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
        data.names = phy$tip.label
    }
    else data.names = rownames(mf)
    order <- match(data.names, phy$tip.label)
    if (sum(is.na(order)) > 0) {
        warning("Data names do not match with the tip labels.\n")
        rownames(mf) <- data.names
    }
    else {
        tmp <- mf
        rownames(mf) <- phy$tip.label
        mf[order, ] <- tmp[1:nrow(tmp), ]
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
    Vphy <- vcv(phy)
    Vphy <- Vphy/max(Vphy)
    Vphy/exp(determinant(Vphy)$modulus[1]/n)
    if (!is.null(B.init) & length(B.init) != p) {
        warning("B.init not correct length, so computed B.init using glm()")
    }
    if (is.null(B.init) | (!is.null(B.init) & length(B.init) != 
        p)) {
        B.init <- t(matrix(glm(formula = formula, data = data, 
            family = "binomial")$coefficients, ncol = p))
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
    while (((t(est.s2 - oldest.s2) %*% (est.s2 - oldest.s2) > 
        tol.pql^2) | (t(est.B - oldest.B) %*% (est.B - oldest.B)/length(B) > 
        tol.pql^2)) & (iteration <= maxit.pql)) {
        iteration <- iteration + 1
        oldest.s2 <- est.s2
        oldest.B <- est.B
        est.B.m <- B
        oldest.B.m <- matrix(10^6, nrow = length(est.B))
        iteration.m <- 0
        while ((t(est.B.m - oldest.B.m) %*% (est.B.m - oldest.B.m)/length(B) > 
            tol.pql^2) & (iteration.m <= maxit.pql)) {
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
        opt <- optim(fn = pglmm.reml, par = s2, tinvW = invW, 
            tH = H, tVphy = Vphy, tX = X, method = "BFGS", control = list(factr = 1e+12, 
                maxit = maxit.reml))
        s2 <- abs(opt$par)
        C <- s2 * Vphy
        est.s2 <- s2
        est.B <- B
    }
    convergeflag <- "converged"
    if (iteration >= maxit.pql | rcondflag >= 3) {
        convergeflag <- "Did not converge; try increasing maxit.pql or starting with B.init values of .001"
    }
    converge.test.s2 <- (t(est.s2 - oldest.s2) %*% (est.s2 - 
        oldest.s2))^0.5
    converge.test.B <- (t(est.B - oldest.B) %*% (est.B - oldest.B))^0.5/length(est.B)
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
    lnlike.cond.reml <- -0.5 * (n - p) * log(2 * pi) + 0.5 * 
        determinant(t(X) %*% X)$modulus[1] - 0.5 * LL
    LL0 <- pglmm.reml(par = 0, tinvW = invW, tH = H, tVphy = Vphy, 
        tX = X)
    lnlike.cond.reml0 <- -0.5 * (n - p) * log(2 * pi) + 0.5 * 
        determinant(t(X) %*% X)$modulus[1] - 0.5 * LL0
    P.H0.s2 <- pchisq(2 * (lnlike.cond.reml - lnlike.cond.reml0), 
        df = 1, lower.tail = F)/2
    results <- list(formula = formula, B = B, B.se = B.se, B.cov = B.cov, 
        B.zscore = B.zscore, B.pvalue = B.pvalue, s2 = s2, P.H0.s2 = P.H0.s2, 
        mu = mu, b = b, B.init = B.init, X = X, y = y, phy = phy, data = data, H = H, VCV = Vphy, 
        V = V, convergeflag = convergeflag, iteration = iteration, 
        converge.test.s2 = converge.test.s2, converge.test.B = converge.test.B, 
        rcondflag = rcondflag)
    class(results) <- "binaryPGLMM"
    results
}

########################################################
########################################################
# version of phyloglm that uses Nealder-Mead optimization
########################################################
########################################################

phyloglm <- function (formula, data = list(), phy, method = c("logistic_MPLE", 
    "logistic_IG10", "poisson_GEE"), btol = 10, log.alpha.bound = 4, 
    start.beta = NULL, start.alpha = NULL, boot = 0, full.matrix = TRUE, opt.method = c("Nelder-Mead","L-BFGS-B","BFGS","SANN")) 
{
    if (!inherits(phy, "phylo")) 
        stop("object \"phy\" is not of class \"phylo\".")
    if (is.null(phy$edge.length)) 
        stop("the tree has no branch lengths.")
    if (is.null(phy$tip.label)) 
        stop("the tree has no tip labels.")
    method = match.arg(method)
    phy = reorder(phy, "pruningwise")
    original.edge.length = phy$edge.length
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- n + 1L
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    mf = model.frame(formula = formula, data = data)
    if (nrow(mf) != length(phy$tip.label)) 
        stop("the number of rows in the data does not match the number of tips in the tree.")
    if (is.null(rownames(mf))) 
        warning("the data has no names, order assumed to be the same as tip labels in the tree.\n")
    else {
        ordr = match(phy$tip.label, rownames(mf))
        if (sum(is.na(ordr)) > 0) 
            stop("the row names in the data do not match the tip labels in the tree.\n")
        mf = mf[ordr, , drop = F]
    }
    X = model.matrix(attr(mf, "terms"), data = mf)
    y = model.response(mf)
    dk = ncol(X)
    dis = pruningwise.distFromRoot(phy)
    if (method %in% c("logistic_MPLE", "logistic_IG10")) {
        if (sum(!(y %in% c(0, 1)))) 
            stop("The model by Ives and Garland requires a binary response (dependent variable).")
        if (var(y) == 0) 
            stop("the response (dependent variable) is always 0 or always 1.")
        btouch = 0
        proposedBetaSD = 0.05
        D = max(dis[1:n]) - dis[1:n]
        D = D - mean(D)
        externalEdge = (des <= n)
        phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + 
            D[des[externalEdge]]
        times <- pruningwise.branching.times(phy)
        names(times) <- (n + 1):(n + phy$Nnode)
        Tmax <- max(times)
        intern = which(phy$edge[, 2] > n)
        lok = rep(-1, N)
        lok[intern] = des[intern] - n
    }
    if (method == "poisson_GEE") {
        if ((!isTRUE(all(y == floor(y))))) 
            stop("The Poisson regression requires an integer response (dependent variable).")
        if (sum(y < 0)) 
            stop("The Poisson regression requires a positive response (dependent variable).")
    }
    transf.branch.lengths_poisson_GEE <- function(beta) {
        if (dk > 1) 
            g = X %*% beta
        else g = rep(1, n) * beta
        mu = as.vector(exp(g))
        root.edge = 0
        diag = sqrt(mu/dis[1:n])
        edge.length = phy$edge.length
        return(list(edge.length, root.edge, diag))
    }
    transf.branch.lengths <- function(B, lL) {
        if (dk > 1) 
            g = X %*% B
        else g = rep(1, n) * B
        mu = as.vector(1/(1 + exp(-g)))
        p = mean(mu)
        alpha = 1/exp(lL)
        edge.length = numeric(N)
        distFromRoot <- exp(-2 * alpha * times)
        tmp = .C("transbranchlengths_IvesGarland2010", as.integer(N), 
            as.integer(des), as.integer(anc - n), as.integer(lok), 
            as.double(distFromRoot), as.integer(externalEdge), 
            as.double(mu), as.double(p), as.double(alpha), as.double(D), 
            el = as.double(1:N), di = as.double(1:n))
        edge.length = tmp$el
        diag = tmp$di
        root.edge = min(distFromRoot)
        if (any(is.nan(edge.length))) 
            stop("edge.length[i] is NaN. Please reduce btol and/or log.alpha.bound.")
        return(list(edge.length, root.edge, diag))
    }
    three.point.compute <- function(trans, y, X) {
        ole = 4 + 2 * dk + dk * dk
        tmp = .C("threepoint", as.integer(N), as.integer(n), 
            as.integer(phy$Nnode), as.integer(1), as.integer(dk), 
            as.integer(ROOT), as.double(trans[[2]]), as.double(trans[[1]]), 
            as.integer(des), as.integer(anc), as.double(as.vector(y)), 
            as.double(as.vector(X)), result = double(ole))$result
        return(list(vec11 = tmp[2], y1 = tmp[3], yy = tmp[4], 
            X1 = tmp[5:(4 + dk)], XX = matrix(tmp[(5 + dk):(ole - 
                dk)], dk, dk), Xy = tmp[(ole - dk + 1):ole], 
            logd = tmp[1]))
    }
    plogregfunct <- function(startB, startlL) {
        convergeflag = 0
        clL = startlL
        cB = startB
        diflL = 100
        difB = 100
        counter = 0
        ttozero = 10^6
        optss <- list(reltol = .Machine$double.eps^0.5, maxit = 1e+05, 
            parscale = 1)
        while (((diflL > 10^-6) | (difB > 10^-6) | (ttozero > 
            10^-1)) & (counter < 20)) {
            counter = counter + 1
            oldlL = clL
            oldB = cB
            olddiflL = diflL
            olddifB = difB
            opt <- optim(par = clL, fn = function(par) {
                plogreglLfunct(cB, par)
            }, method = opt.method)
            clL = as.numeric(opt$par)
            diflL = (clL - oldlL)^2
            if (counter >= 10) 
                clL = (clL + oldlL)/2
            opt <- optim(par = cB, fn = function(par) {
                plogregBfunct(par, clL)
            }, method = opt.method, control = list(factr = 1e+12))
            cB = as.vector(opt$par)
            ttozero = as.numeric(opt$value)
            if (ttozero > 10^-2) {
                Btemp = rnorm(dk, startB, proposedBetaSD * pmax(abs(startB), 
                  rep(0.1, dk)))
                opt <- optim(par = Btemp, fn = function(par) {
                  plogregBfunct(par, clL)
                }, method = opt.method, control = list(factr = 1e+12))
                Btemp = as.vector(opt$par)
                newttozero = as.numeric(opt$value)
                if (newttozero < ttozero) {
                  cB = Btemp
                  ttozero = newttozero
                }
            }
            difB = sum((cB - oldB) * (cB - oldB))
            if (counter >= 10) 
                cB = (cB + oldB)/2
        }
        if (counter >= 19) 
            if ((max(abs(c(oldlL - clL, oldB - cB))) > 0.1) | 
                (ttozero > 0.5)) 
                convergeflag = 1
        return(list(B = cB, lL = clL, convergeflag = convergeflag))
    }
    plogregBfunct <- function(B, lL) {
        if (dk > 1) 
            g = X %*% B
        else g = rep(1, n) * B
        if (any(abs(g) >= btol)) {
            btouch <<- 1
            return(1e+06)
        }
        mu = as.vector(1/(1 + exp(-g)))
        temp = transf.branch.lengths(B, lL)
        dia = temp[[3]]
        comp = three.point.compute(temp[1:2], (y - mu)/dia, mu * 
            (1 - mu) * X/dia)
        logdetC = comp$logd + 2 * sum(log(dia)) - sum(log(mu * 
            (1 - mu)))
        if (logdetC < -100 * log(10)) 
            return(1e+06)
        Z = comp$Xy
        if (dk == 1) 
            FirthC = (1 - 2 * mu)/2
        else {
            Dx = 0.1
            infoM = comp$XX
            invInfoM = solve(infoM)
            FirthC = rep(NA, dk)
            for (i in 1:dk) {
                dB = B
                dB[i] = dB[i] + Dx
                g = X %*% dB
                if (any(abs(g) >= btol)) 
                  return(1e+06)
                mu = as.vector(1/(1 + exp(-g)))
                ttemp = transf.branch.lengths(dB, lL)
                tdiag = ttemp[[3]]
                tcomp = three.point.compute(ttemp[1:2], (y - 
                  mu)/tdiag, mu * (1 - mu) * X/tdiag)
                dinfoMp = tcomp$XX
                dB = B
                dB[i] = dB[i] - Dx
                g = X %*% dB
                if (any(abs(g) >= btol)) 
                  return(1e+06)
                mu = as.vector(1/(1 + exp(-g)))
                ttemp = transf.branch.lengths(dB, lL)
                tdiag = ttemp[[3]]
                tcomp = three.point.compute(ttemp[1:2], (y - 
                  mu)/tdiag, mu * (1 - mu) * X/tdiag)
                dinfoMm = tcomp$XX
                DinfoM = (dinfoMp - dinfoMm)/Dx/2
                FirthC[i] = sum(diag(invInfoM %*% DinfoM))/2
            }
        }
        tozero = Z + FirthC
        return(sum(tozero^2))
    }
    plogreglLfunct <- function(B, lL) {
        g = X %*% B
        mu = as.vector(1/(1 + exp(-g)))
        if (abs(lL - log(Tmax)) >= log.alpha.bound) 
            return(1e+10)
        temp = transf.branch.lengths(B, lL)
        dia = temp[[3]]
        comp = three.point.compute(temp[1:2], (y - mu)/dia, mu * 
            (1 - mu) * X/dia)
        LL = (comp$logd + 2 * sum(log(dia)) + comp$yy)/2
        if (!is.finite(LL)) 
            LL = 1e+10
        return(LL)
    }
    plogregBSEfunct <- function(B, lL) {
        g = X %*% B
        mu = as.vector(1/(1 + exp(-g)))
        temp = transf.branch.lengths(B, lL)
        dia = temp[[3]]
        comp = three.point.compute(temp[1:2], (y - mu)/dia, mu * 
            (1 - mu) * X/dia)
        infoM = comp$XX
        covBSE = solve(infoM)
        BSE = sqrt(diag(covBSE))
        return(list(BSE = BSE, covBSE = covBSE, info = infoM))
    }
    npllh <- function(par) {
        if (abs(par[dk + 1] - log(Tmax)) >= log.alpha.bound) 
            return(1e+10)
        g = X %*% par[1:dk]
        if (any(abs(g) >= btol)) {
            btouch <<- 1
            return(1e+10)
        }
        mu = as.vector(1/(1 + exp(-g)))
        temp = transf.branch.lengths(par[1:dk], par[dk + 1])
        dia = temp[[3]]
        comp = three.point.compute(temp[1:2], numeric(n), mu * 
            (1 - mu) * X/dia)
        infoM = comp$XX
        llk <- .C("logistreglikelihood", as.integer(N), as.integer(n), 
            as.integer(phy$Nnode), as.integer(ROOT), as.double(original.edge.length), 
            as.integer(des), as.integer(anc), as.integer(as.vector(y)), 
            as.double(as.vector(mu)), as.integer(dk), as.double(exp(-par[dk + 
                1])), loglik = double(1))$loglik
        if (dk == 1) 
            pllik = llk + log(abs(infoM))/2
        else pllik = llk + log(det(infoM))/2
        -pllik
    }
    llh <- function(mu, alpha) {
        .C("logistreglikelihood", as.integer(N), as.integer(n), 
            as.integer(phy$Nnode), as.integer(ROOT), as.double(original.edge.length), 
            as.integer(des), as.integer(anc), as.integer(as.vector(y)), 
            as.double(as.vector(mu)), as.integer(dk), as.double(alpha), 
            loglik = double(1))$loglik
    }
    iterate_beta <- function(beta) {
        difbeta = 1
        maxint = 10000
        count = 0
        curbeta = beta
        while ((difbeta > 1e-10) && (count < maxint)) {
            mu = as.vector(exp(X %*% curbeta))
            temp = transf.branch.lengths_poisson_GEE(curbeta)
            dia = temp[[3]]
            if (sum(which(mu == 0)) > 0) 
                break
            comp = three.point.compute(temp[1:2], (y - mu)/dia, 
                mu * X/dia)
            invI = solve(comp$XX)
            newbeta = curbeta + invI %*% comp$Xy
            count = count + 1
            difbeta = sum(abs(newbeta - curbeta))
            curbeta = newbeta
        }
        mu = as.vector(exp(X %*% curbeta))
        r = (y - mu)/sqrt(mu)
        phi = sum(r^2)/(n - dk)
        covBSE = phi * invI
        BSE = sqrt(diag(covBSE))
        if (difbeta > 1e-10) 
            convergeflag = 1
        else convergeflag = 0
        return(list(beta = as.vector(curbeta), BSE = BSE, covBSE = covBSE, 
            phi = phi, convergeflag = convergeflag))
    }
    if (is.null(start.beta)) {
        if (method %in% c("logistic_MPLE", "logistic_IG10")) {
            fit = glm(y ~ X - 1, family = binomial)
            startB = fit$coefficients
            if (any(abs(X %*% startB) >= btol)) {
                warning("The estimated coefficients in the absence of phylogenetic signal lead\n  to some linear predictors beyond 'btol'. Increase btol?\n  Starting from beta=0 other than intercept.")
                startB = numeric(dk)
                iint = match("(Intercept)", colnames(X))
                if (!is.na(iint)) 
                  startB[iint] = log(sum(y == 1)/sum(y == 0))
                if (any(abs(X %*% startB) >= btol)) 
                  startB[iint] = 0
            }
        }
        if (method == "poisson_GEE") {
            fit = glm(y ~ X - 1, family = poisson)
            start.beta = fit$coefficients
        }
    }
    else {
        if (length(start.beta) != dk) 
            stop(paste("start.beta shoudl be of length", dk))
        if (method %in% c("logistic_MPLE", "logistic_IG10")) {
            startB = as.vector(start.beta)
            if (any(abs(X %*% startB) >= btol)) 
                stop("With these starting beta values, some linear predictors are beyond 'btol'.\n  Increase btol or choose new starting values for beta.")
        }
    }
    if (method %in% c("logistic_MPLE", "logistic_IG10")) {
        if (is.null(start.alpha)) 
            startlL = log(Tmax)
        else {
            if (length(start.alpha) != 1) 
                stop("start.alpha should be a single positive value")
            if (start.alpha <= 0) 
                stop("start.alpha should be a positive value")
            startlL = -log(start.alpha)
            if (abs(startlL - log(Tmax)) >= log.alpha.bound) {
                tmp = "start.alpha is outside the bounds, which are\n  exp(+/-log.alpha.bound)/Tmax: "
                tmp = paste(tmp, signif(exp(-log.alpha.bound)/Tmax, 
                  3), ",", signif(exp(log.alpha.bound)/Tmax, 
                  3), " (Tmax=", Tmax, ").", "\n  Change start.alpha or increase log.alpha.bound.", 
                  sep = "")
                stop(tmp)
            }
        }
    }
    if (method %in% c("logistic_MPLE", "logistic_IG10")) {
        if (method == "logistic_IG10") {
            plogreg = plogregfunct(startB, startlL)
            lL = plogreg$lL
            B = plogreg$B
            convergeflag = plogreg$convergeflag
        }
        if (method == "logistic_MPLE") {
            opt <- optim(par = c(startB, startlL), fn = npllh, 
                method = opt.method, control = list(factr = 1e+12))
            B = opt$par[1:dk]
            lL = opt$par[dk + 1]
            convergeflag = opt$convergence
        }
        if ((lL - log(Tmax) + 0.02) > log.alpha.bound) {
            warn = paste("the estimate of 'alpha' (", 1/exp(lL), 
                ") reached the lower bound (", 1/Tmax/exp(log.alpha.bound), 
                ").\n This may reflect a flat likelihood at low alpha values near 0,\n", 
                " meaning that the phylogenetic correlation is estimated to be maximal\n", 
                " under the model in Ives and Garland (2010).", 
                sep = "")
            warning(warn)
        }
        if ((lL - log(Tmax) - 0.02) < -log.alpha.bound) {
            warn = paste("the estimate of 'alpha' (", 1/exp(lL), 
                ") reached the upper bound (", exp(log.alpha.bound)/Tmax, 
                ").\n This may simply reflect a flat likelihood at large alpha values,\n", 
                " meaning that the phylogenetic correlation is estimated to be negligible.", 
                sep = "")
            warning(warn)
        }
        if (btouch == 1) 
            warning("the boundary of the linear predictor has been reached during the optimization procedure.\nYou can increase this bound by increasing 'btol'.")
        plogregBSE = plogregBSEfunct(B, lL)
        results <- list(coefficients = B, alpha = 1/exp(lL), 
            sd = plogregBSE$BSE, vcov = plogregBSE$covBSE, convergence = convergeflag)
    }
    if (method == "poisson_GEE") {
        res = iterate_beta(as.vector(start.beta))
        results <- list(coefficients = res$beta, scale = res$phi, 
            sd = res$BSE, vcov = res$covBSE, convergence = res$convergeflag)
    }
    if (results$converge) 
        warning("phyloglm failed to converge.\n")
    names(results$coefficients) = colnames(X)
    colnames(results$vcov) = colnames(X)
    rownames(results$vcov) = colnames(X)
    results$linear.predictors = as.vector(X %*% results$coefficients)
    names(results$linear.predictors) = names(y)
    if (method %in% c("logistic_MPLE", "logistic_IG10")) {
        if (max(abs(results$linear.predictors)) + 0.01 > btol) 
            warning("the linear predictor reaches its bound for one (or more) tip.")
        results$fitted.values = as.vector(1/(1 + exp(-results$linear.predictors)))
        results$mean.tip.height = Tmax
        results$logLik = llh(results$fitted.values, results$alpha)
        results$penlogLik = results$logLik + log(det(as.matrix(plogregBSE$info)))/2
        results$aic = -2 * results$logLik + 2 * (dk + 1)
    }
    if (method == "poisson_GEE") {
        results$fitted.values = as.vector(exp(-results$linear.predictors))
        results$logLik = NA
        results$penlogLik = NA
        results$aic = NA
    }
    names(results$fitted.values) = names(y)
    results$residuals = y - results$fitted.values
    results$y = y
    results$n = n
    results$d = dk
    results$formula = formula
    results$call = match.call()
    results$method = method
    results$X = X
    results$boot = boot
    if ((boot > 0) && (method %in% c("logistic_MPLE", "logistic_IG10"))) {
        options(warn = -1)
        bootobject <- rbinTrait(n = boot, phy = phy, beta = results$coefficients, 
            alpha = results$alpha, X = X, model = "LogReg")
        ncoeff = length(results$coefficients)
        bootmatrix <- matrix(NA, boot, ncoeff + 1)
        colnames(bootmatrix) <- c(names(results$coefficients), 
            "alpha")
        for (i in 1:boot) {
            y = bootobject[, i]
            if (method == "logistic_IG10") {
                bootfit <- try(plogregfunct(startB, startlL), 
                  silent = TRUE)
                if (!inherits(bootfit, "try-error")) {
                  bootmatrix[i, 1:ncoeff] <- bootfit$B
                  bootmatrix[i, ncoeff + 1] <- 1/exp(bootfit$lL)
                }
            }
            if (method == "logistic_MPLE") {
                bootfit <- try(optim(par = c(startB, startlL), 
                  fn = npllh, method = opt.method, control = list(factr = 1e+12)), 
                  silent = TRUE)
                if (!inherits(bootfit, "try-error")) {
                  bootmatrix[i, 1:ncoeff] <- bootfit$par[1:dk]
                  bootmatrix[i, ncoeff + 1] <- 1/exp(bootfit$par[dk + 
                    1])
                }
            }
        }
        ind.na <- which(is.na(bootmatrix[, 1]))
        if (length(ind.na) > 0) {
            bootmatrix <- bootmatrix[-ind.na, ]
            numOnes <- range(apply(bootobject[, ind.na], 2, sum))
        }
        bootmean <- apply(bootmatrix, 2, mean)
        bootsd <- apply(bootmatrix, 2, sd)
        bootconfint95 <- apply(bootmatrix, 2, quantile, probs = c(0.025, 
            0.975))
        bootmeanAlog <- mean(log(bootmatrix[, ncoeff + 1]))
        bootsdAlog <- sd(log(bootmatrix[, ncoeff + 1]))
        results$bootmean = bootmean
        results$bootsd = bootsd
        results$bootconfint95 = bootconfint95
        results$bootmeanAlog = bootmeanAlog
        results$bootsdAlog = bootsdAlog
        results$bootnumFailed = length(ind.na)
        if (full.matrix) 
            results$bootstrap = bootmatrix
        options(warn = 0)
    }
    class(results) = "phyloglm"
    results
}
