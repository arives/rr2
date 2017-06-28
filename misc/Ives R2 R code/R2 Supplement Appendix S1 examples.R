# Illustrations of R2.ls, R2.lr, and R2.ce for LMM, PGLS, GLMMs for binary (binomial) data, and PLOG. Data are simulated using the same models used in the main text.

# R2.ls, R2.lr, and R2.ce functions are called for full model mod and reduced model mod.r as

# R2.ls(mod, mod.r)

# If mod.r is missing, the full R2 is calculated with mod.r as the model with only an intercept.

# mod must be of class lmerMod, glmerMod (from the lme4 package), phylolm and phyloglm (from the phyloglm package), or binaryPGLMM (from the ape package). In addition, mod.r can be class lm and glm.

# The file "R2_Supplement_Appendix_S1_source_code.R" contains the code for R2.ls, R2.lr, and R2.ce; it should be placed in a path that is accessible. IThe file also contains a modified version of phyloglm that performs optimization with Nelder-Mead, and a modified version of binaryPGLMM that includes the dependent variable as an attribute. (I'll ask Emmanuel Paradis to update this in ape when he does the next version.) Because of these modifications, "R2_Supplement_Appendix_S1_source_code.R" should be loaded after ape and phylolm.

# library(phylolm)
# library(ape)
# library(lme4)
library(rr2)
# source("R2_Supplement_Appendix_S1_source_code.R")

#################
# LMM
p1 <- 10
nsample <- 10
n <- p1 * nsample

d <- data.frame(x1=0, x2=0, y=0, u1=rep(1:p1, each=nsample), u2=rep(1:p1, times=nsample))
d$u1 <- as.factor(d$u1)
d$u2 <- as.factor(d$u2)

b1 <- 1
b2 <- -1
sd1 <- 1.5
set.seed(123)
d$x1 <- rnorm(n=n)
d$x2 <- rnorm(n=n)
d$y <- b1 * d$x1 + b2 * d$x2 + rep(rnorm(n=p1, sd=sd1), each=nsample) + rep(rnorm(n=p1, sd=sd1), times=nsample) + rnorm(n=n)

z.f <- lme4::lmer(y ~ x1 + x2 + (1 | u1) + (1 | u2), data=d, REML = F)
z.x <- lme4::lmer(y ~ x1 + (1 | u1) + (1 | u2), data=d, REML = F)
z.v <- lme4::lmer(y ~ 1 + (1 | u2), data=d, REML = F)
z.0 <- lm(y ~ 1, data=d)

R2.ls(z.f, z.x)
R2.ls(z.f, z.v)
R2.ls(z.f)

R2.ce(z.f, z.x)
R2.ce(z.f, z.v)
R2.ce(z.f)

R2.lr(z.f, z.x)
R2.lr(z.f, z.v)
R2.lr(z.f)

# These give the same results
R2.ls(z.f, z.0)
R2.ls(z.f)

# REML is updated to ML for R2.lr, but not R2.ls or R2.ce
z.f <- lme4::lmer(y ~ x1 + x2 + (1 | u1) + (1 | u2), data=d, REML = T)
R2.lr(z.f)
R2.ls(z.f)
R2.ce(z.f)

#################
# PGLS

n <- p1 * nsample
d <- data.frame(x=array(0, dim=n), y=0)

b1 <- 1.5
signal <- 0.7

phy <- ape::compute.brlen(ape::rtree(n=n), method = "Grafen", power = 1)
phy.x <- phy

x <- ape::rTraitCont(phy.x, model = "BM", sigma = 1)

e <- signal^0.5 * ape::rTraitCont(phy, model = "BM", sigma = 1) + (1-signal)^0.5 * rnorm(n=n)
d$x <- x[match(names(e), names(x))]

d$y <- b1 * x + e
rownames(d) <- phy$tip.label

z.x <- phylolm::phylolm(y ~ 1, phy=phy, data=d, model="lambda")
lam.x <- round(z.x$optpar, digits=4)
z.f <- phylolm::phylolm(y ~ x, phy=phy, data=d, model="lambda", starting.value=.98*lam.x+.01)
z.v <- lm(y ~ x, data=d)

R2.ls(z.f, z.x, phy = phy)
R2.ls(z.f, z.v, phy = phy)
R2.ls(z.f, phy = phy)

R2.ce(z.f, z.x, phy = phy)
R2.ce(z.f, z.v, phy = phy)
R2.ce(z.f, phy = phy)

R2.lr(z.f, z.x)
R2.lr(z.f, z.v)
R2.lr(z.f)

#################
# GLMM
library(lme4)

inv.logit <- function(x){
	1/(1 + exp(-x))
}

p1 <- 10
nsample <- 10
n <- p1 * nsample

d <- data.frame(x=0, y=0, u=rep(1:p1, each=nsample))
d$u <- as.factor(d$u)

b1 <- 1
sd1 <- 1.5

d$x <- rnorm(n=n)
prob <- inv.logit(b1 * d$x + rep(rnorm(n=p1, sd=sd1), each=nsample))
d$y <- rbinom(n=n, size=1, prob=prob)

z.f <- lme4::glmer(y ~ x + (1 | u), data=d, family="binomial")
z.x <- lme4::glmer(y ~ 1 + (1 | u), data=d, family="binomial")
z.v <- glm(y ~ x, data=d, family="binomial")

R2.ls(z.f, z.x)
R2.ls(z.f, z.v)
R2.ls(z.f)

system.time(t1 <- R2.ce(z.f, z.x))
system.time(t2 <- R2.ce(z.f, z.v))
system.time(t2.rr2 <- rr2::R2.ce(z.f, z.v))
system.time(t3 <- R2.ce(z.f))
system.time(t3.rr2 <- rr2::R2.ce(z.f))
t1
t2; t2.rr2
t3; t3.rr2

R2.lr(z.f, z.x)
R2.lr(z.f, z.v)
R2.lr(z.f)


#################
# PLOG
set.seed(123)
n <- p1 * nsample
b1 <- 1.5
signal <- 2

phy <- ape::compute.brlen(ape::rtree(n=n), method = "Grafen", power = 1)

# Generate random data
x <- rnorm(n)
d$x <- x

e <- signal * ape::rTraitCont(phy, model = "BM", sigma = 1)
e <- e[match(phy$tip.label, names(e))]

d$y <- rbinom(n=n, size=1, prob=rr2::inv.logit(b1 * d$x + e))
rownames(d) <- phy$tip.label

# R2.ls R2.ce
z.f <- rr2::binaryPGLMM(formula = "y ~ x", data=d, phy=phy)
z.x <- rr2::binaryPGLMM(y ~ 1, data=d, phy=phy)
z.v <- glm(y ~ x, data=d, family="binomial")

R2.ls(z.f, z.x)
rr2::R2.ls(z.f, z.x)
R2.ls(z.f, z.v)
rr2::R2.ls(z.f, z.v)
R2.ls(z.f)
rr2::R2.ls(z.f)


R2.ce(z.f, z.x)
R2.ce(z.f, z.v)
R2.ce(z.f)

# R.lr
z.f <- phyloglm(y ~ x, data=d, start.alpha = 1, phy=phy, opt.method="Nelder-Mead")
z.x <- phyloglm(y ~ 1, data=d, phy=phy, start.alpha=min(20,z.f$alpha), opt.method="Nelder-Mead")
z.v <- glm(y ~ x, data=d, family="binomial")

R2.lr(z.f, z.x)
rr2::R2.lr(z.f, z.x)
R2.lr(z.f, z.v)
rr2::R2.lr(z.f, z.v)
R2.lr(z.f)
rr2::R2.lr(z.f)

